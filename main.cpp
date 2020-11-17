#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <regex>

using namespace std;

typedef struct Pin {
    string name;
    string direct;
    double cap{};
} Pin;

typedef struct CellLibrary {
    string name;
    map<string, Pin> pins;
    vector<double> index_1, index_2;
    vector<vector<double>> cell_rise;
    vector<vector<double>> cell_fall;
    vector<vector<double>> rise_transition;
    vector<vector<double>> fall_transition;
} CellLibrary;

typedef struct Gate {
    string name;
    string type;
    vector<string> input_ports;
    vector<string> input_wires;
    string output_port;
    string output_wire;
    int input_num_in_topological{};
    bool isOutput{};
    list<int> next;
    list<int> prev;
    double critical_cell_delay{};
    double critical_transition{};
} Gate;

class Module {
public:
    string name;
    vector<string> input;
    vector<string> output;
    vector<string> wire;
    vector<Gate> gates;
    vector<int> output_gates;
    vector<int> topological_order;
    vector<double> max_delay;
    vector<int> max_delay_path_from_prev_gate;
    int critical_path_out_index{};

    void topological_sort();

    void delay(map<string, CellLibrary> &cells);

    void output_file();

private:
    void max_delay_path_calculation(int index, double total_output_net_cap, CellLibrary &cell);

    static double
    lookup_table(double &row, double &col, vector<double> &row_index, vector<double> &col_index,
                 vector<vector<double> > &table);
};

void Module::topological_sort() {
    list<int> queue;
    for (int i = 0; i < gates.size(); i++) {
        gates[i].input_num_in_topological = 0;
        gates[i].isOutput = false;

        // Find and flag output gates
        for (auto &j : output) {
            if (std::find(gates[i].input_wires.begin(), gates[i].input_wires.end(), j) != gates[i].input_wires.end()) {
                output_gates.push_back(i);
                gates[i].isOutput = true;
                break;
            }
        }
        // Find and flag input gates
        for (auto &j : input) {
            if (std::find(gates[i].input_wires.begin(), gates[i].input_wires.end(), j) != gates[i].input_wires.end()) {
                gates[i].input_num_in_topological++;
            }
        }

        if (gates[i].input_num_in_topological == gates[i].input_wires.size()) {
            // Doesn't have input from other gates
            queue.push_back(i);
        }
    }


    while (!queue.empty()) {
        int index = queue.front();
        queue.pop_front();
        topological_order.push_back(index);

        string out = gates[index].output_wire;
        for (int i = 0; i < gates.size(); i++) {
            for (int j = 0; j < gates[i].input_wires.size(); j++) {
                if (gates[i].input_wires[j] == out) {
                    gates[i].input_num_in_topological++;
                    gates[i].prev.push_back(index);
                    gates[index].next.push_back(i);
                    if (gates[i].input_num_in_topological == gates[i].input_wires.size()) {
                        queue.push_back(i);
                    }
                    break;
                }
            }
        }
    }
}

double Module::lookup_table(double &row, double &col, vector<double> &row_index, vector<double> &col_index,
                            vector<vector<double> > &table) {
    int row_end_index = (int) row_index.size() - 1;
    int col_end_index = (int) col_index.size() - 1;
    int index_row1 = row_end_index;
    int index_row2 = 0;
    int index_col1 = col_end_index;
    int index_col2 = 0;

    for (int i = row_end_index; i >= 0; i--) if (row_index[i] > row) index_row1--;
    for (int i = 0; i <= row_end_index; i++) if (row_index[i] <= row) index_row2++;
    for (int i = col_end_index; i >= 0; i--) if (col_index[i] > col) index_col1--;
    for (int i = 0; i <= col_end_index; i++) if (col_index[i] <= col) index_col2++;


    if (index_row1 < 0) {
        index_row1++;
        index_row2++;
    }
    if (index_row2 > row_end_index) {
        index_row1--;
        index_row2--;
    }
    if (index_col1 < 0) {
        index_col1++;
        index_col2++;
    }
    if (index_col2 > col_end_index) {
        index_col1--;
        index_col2--;
    }

    double col2_col = col_index[index_col2] - col;
    double col2_col1 = col_index[index_col2] - col_index[index_col1];
    double col_col1 = col - col_index[index_col1];
    double row2_row = row_index[index_row2] - row;
    double row2_row1 = row_index[index_row2] - row_index[index_row1];
    double row_row1 = row - row_index[index_row1];
    return col2_col * row2_row / col2_col1 / row2_row1 * table[index_row1][index_col1]
           + col_col1 * row2_row / col2_col1 / row2_row1 * table[index_row1][index_col2]
           + col2_col * row_row1 / col2_col1 / row2_row1 * table[index_row2][index_col1]
           + col_col1 * row_row1 / col2_col1 / row2_row1 * table[index_row2][index_col2];
}

void Module::max_delay_path_calculation(int cur_gate_index, double total_output_net_cap, CellLibrary &cell) {
    int max_path_prev_gate_index = -1;
    double prev_max_delay = 0;
    double prev_max_transition = 0;
    if (!gates[cur_gate_index].prev.empty()) {
        prev_max_delay = max_delay[gates[cur_gate_index].prev.front()];
        prev_max_transition = gates[gates[cur_gate_index].prev.front()].critical_transition;
        max_delay_path_from_prev_gate[cur_gate_index] = gates[cur_gate_index].prev.front();

        for (auto lp = gates[cur_gate_index].prev.begin();
             lp != gates[cur_gate_index].prev.end(); lp++) {
            if (gates[*lp].critical_transition > prev_max_transition) {
                prev_max_transition = gates[*lp].critical_transition;
            }
            if (max_delay[*lp] > prev_max_delay) {
                prev_max_delay = max_delay[*lp];
                max_delay_path_from_prev_gate[cur_gate_index] = *lp;
            }
        }
    }
    max_path_prev_gate_index = max_delay_path_from_prev_gate[cur_gate_index];

    double max_input_transition = prev_max_transition;

    double cell_rise = lookup_table(max_input_transition, total_output_net_cap, cell.index_2, cell.index_1,
                                    cell.cell_rise);
    double cell_fall = lookup_table(max_input_transition, total_output_net_cap, cell.index_2, cell.index_1,
                                    cell.cell_fall);
    gates[cur_gate_index].critical_cell_delay = (cell_rise > cell_fall) ? cell_rise : cell_fall;

    double rise_transition = lookup_table(max_input_transition, total_output_net_cap, cell.index_2, cell.index_1,
                                          cell.rise_transition);
    double fall_transition = lookup_table(max_input_transition, total_output_net_cap, cell.index_2, cell.index_1,
                                          cell.fall_transition);
    gates[cur_gate_index].critical_transition = (rise_transition + fall_transition) / 2.0;

    max_delay[cur_gate_index] = prev_max_delay + gates[cur_gate_index].critical_cell_delay;
}

// TODO: Fix stack error
void Module::delay(map<string, CellLibrary> &cells) {
    for (int i = 0; i < gates.size(); i++) {
        max_delay.push_back(0);
        max_delay_path_from_prev_gate.push_back(-1);
    }

    for (int cur_gate_index : topological_order) {
        double total_output_net_cap = (gates[cur_gate_index].isOutput) ? 0.03 : 0.0; // Primary output loading is 0.03
        for (auto lp = gates[cur_gate_index].next.begin();
             lp != gates[cur_gate_index].next.end(); lp++) {
            for (int k = 0; k < gates[*lp].input_wires.size() - 1; k++) {
                if (gates[cur_gate_index].input_wires.back() == gates[*lp].input_wires[k]) {
                    total_output_net_cap += cells[gates[*lp].type].pins[gates[*lp].input_ports[k]].cap;
                    break;
                }
            }
        }

        max_delay_path_calculation(cur_gate_index, total_output_net_cap, cells[gates[cur_gate_index].type]);
    }

    critical_path_out_index = output_gates.front();
    for (int i = 1; i < output_gates.size(); i++) {
        if (max_delay[output_gates[i]] > max_delay[critical_path_out_index])
            critical_path_out_index = output_gates[i];
    }
}

void Module::output_file() {
    fstream fout;
    fout.open("timing.txt", ios::out);
    for (auto &gate : gates) {
        fout << gate.name << " "
             << gate.critical_cell_delay << " "
             << gate.critical_transition << endl;
    }
    fout.close();


    fout.open("critical_delay_path.txt", ios::out);
    fout << "Longest delay = " << max_delay[critical_path_out_index] << ", the path is:" << endl;
    vector<string> max_delay_path;
    int index = critical_path_out_index;
    while (max_delay_path_from_prev_gate[index] != -1) {
        max_delay_path.push_back(gates[index].input_wires.back());
        index = max_delay_path_from_prev_gate[index];
    }
    max_delay_path.push_back(gates[index].input_wires.back());
    max_delay_path.push_back(gates[index].input_wires.front());
    for (int i = max_delay_path.size() - 1; i > 0; --i)
        fout << max_delay_path[i] << " -> ";
    fout << max_delay_path.front() << endl;
    fout.close();
}

enum line {
    module,
    input,
    output,
    wire,
    others
};

class Parser {
public:
    static void parse_cell_lib(char *file, map<string, CellLibrary> &cells);

    static void parse_module(char *file, Module &m, const map<string, CellLibrary> &cells);

private:
    static void parse_lu_table(fstream &fin, vector<double> &i1, vector<double> &i2);

    static void parse_cell(fstream &fin, CellLibrary &c);

    static void parse_pin(fstream &fin, Pin &p);

    static void parse_cell_rise(fstream &fin, CellLibrary &c);

    static void parse_cell_fall(fstream &fin, CellLibrary &c);

    static void parse_rise_transition(fstream &fin, CellLibrary &c);

    static void parse_fall_transition(fstream &fin, CellLibrary &c);

    static vector<string> tokenize(string &raw_str);

    static string remove_comment_and_space(string &original_str);
};

void Parser::parse_cell_lib(char *file, map<string, CellLibrary> &cells) {
    fstream fin;
    fin.open(file, ios::in);
    if (!fin) {
        cerr << "Fail to open: " << file << endl;
        exit(1);
    }

    string in, temp;
    int num;
    stringstream ss;
    vector<double> index_1, index_2;

    while (getline(fin, in)) {
        in = remove_comment_and_space(in);
        ss << in;
        while (getline(ss, temp, '(')) {
            if (temp == "lu_table_template") {
                parse_lu_table(fin, index_1, index_2);
            } else if (temp == "cell") {
                CellLibrary c;
                c.index_1 = index_1;
                c.index_2 = index_2;
                getline(ss, c.name, ')');
                parse_cell(fin, c);
                cells[c.name] = c;
            }
        }

        ss.str("");
        ss.clear();
    }

    fin.close();
}

vector<string> Parser::tokenize(string &raw_str) {
    // Remove comments
    string cleaned_str;
    const std::regex comment(R"(\/\/[^\n\r]*\n|\/\*[\s\S]*\*\/)");    // Verilog comment format
    std::regex_replace(std::back_inserter(cleaned_str), raw_str.begin(), raw_str.end(), comment, "");

    // Parse code
    const std::regex re(R"([\s\t\n\(\)|,]+)");
    std::sregex_token_iterator it{cleaned_str.begin(), cleaned_str.end(), re, -1};
    std::vector<std::string> tokenized{it, {}};
    tokenized.erase(
            std::remove_if(tokenized.begin(),
                           tokenized.end(),
                           [](std::string const &s) {
                               return s.empty();
                           }),
            tokenized.end());

    return tokenized;
}

string Parser::remove_comment_and_space(string &original_str) {
    string cleaned_str;
    const std::regex reg(R"(\/\/[^\n\r]*\n|\/\*[\s\S]*\*\/|[\s\t\n]*)");    // Verilog comment format
    std::regex_replace(std::back_inserter(cleaned_str), original_str.begin(), original_str.end(), reg, "");
    return cleaned_str;
}

void Parser::parse_lu_table(fstream &fin, vector<double> &i1, vector<double> &i2) {
    string in, temp;
    stringstream ss;

    const std::regex idx1_rgx("index_1\\(\"([^\\;\\)]+)\"\\)");
    const std::regex idx2_rgx("index_2\\(\"([^\\;\\)]+)\"\\)");
    // Extract lu_table_template content
    while (getline(fin, temp, ';')) {
        temp = remove_comment_and_space(temp);
        std::smatch results;
        if (std::regex_search(temp, results, idx1_rgx)) {
            temp = (results.begin() + 1)->str();    // temp is the value between ""
            ss.clear();
            ss << temp;
            while (getline(ss, temp, ',')) {
                i1.push_back(strtod(temp.c_str(), nullptr));
            }
        }
        if (std::regex_search(temp, results, idx2_rgx)) {
            temp = (results.begin() + 1)->str();
            ss.clear();
            ss << temp;
            while (getline(ss, temp, ',')) {
                i2.push_back(strtod(temp.c_str(), nullptr));
            }
            break;
        }
    }
}

void Parser::parse_cell(fstream &fin, CellLibrary &c) {
    string in, temp;
    stringstream ss1, ss2;
    int cnt = 0;

    while (cnt < 4) {   // 4 timing prooerity each cell
        getline(fin, in);
        ss2 << in;
        ss2 >> temp;
        ss1 << temp;
        getline(ss1, temp, '(');
        temp = remove_comment_and_space(temp);
        if (temp == "pin") {
            Pin p;
            getline(ss1, p.name, ')');
            parse_pin(fin, p);
            c.pins[p.name] = p;
        } else if (temp == "cell_rise") {
            parse_cell_rise(fin, c);
            cnt++;
        } else if (temp == "cell_fall") {
            parse_cell_fall(fin, c);
            cnt++;
        } else if (temp == "rise_transition") {
            parse_rise_transition(fin, c);
            cnt++;
        } else if (temp == "fall_transition") {
            parse_fall_transition(fin, c);
            cnt++;
        }

        ss1.str("");
        ss1.clear();
        ss2.str("");
        ss2.clear();
    }
}

void Parser::parse_pin(fstream &fin, Pin &p) {
    string in, temp;
    const std::regex input_rgx("direction:input");
    const std::regex output_rgx("direction:output");
    const std::regex capacitance_rgx("capacitance:([^;]+)");
    std::smatch results;

    int cnt = 0;
    while (cnt < 2) {
        getline(fin, temp);
        temp = remove_comment_and_space(temp);
        if (std::regex_search(temp, results, input_rgx)) {
            p.direct = "input";
            ++cnt;
        } else if (std::regex_search(temp, results, output_rgx)) {
            p.direct = "output";
            ++cnt;
        } else if (std::regex_search(temp, results, capacitance_rgx)) {
            temp = (results.begin() + 1)->str();    // temp is the value between ""
            p.cap = strtod(temp.c_str(), nullptr);
            ++cnt;
        }
    }
}

void Parser::parse_cell_rise(fstream &fin, CellLibrary &c) {
    string in, temp;
    int first, last;
    stringstream ss;

    for (int i = 0; i < c.index_1.size(); i++) {
        getline(fin, temp);
        first = temp.find('"');
        last = temp.find_last_of('"');
        in = temp.substr(first + 1, last - first - 1);
        ss << in;

        vector<double> value;
        while (getline(ss, temp, ',')) {
            value.push_back(atof(temp.c_str()));
        }
        c.cell_rise.push_back(value);

        ss.str("");
        ss.clear();
    }
}

void Parser::parse_cell_fall(fstream &fin, CellLibrary &c) {
    string in, temp;
    int first, last;
    stringstream ss;

    for (int i = 0; i < c.index_1.size(); i++) {
        getline(fin, temp);
        first = temp.find('"');
        last = temp.find_last_of('"');
        in = temp.substr(first + 1, last - first - 1);
        ss << in;

        vector<double> value;
        while (getline(ss, temp, ',')) {
            value.push_back(atof(temp.c_str()));
        }
        c.cell_fall.push_back(value);

        ss.str("");
        ss.clear();
    }
}

void Parser::parse_rise_transition(fstream &fin, CellLibrary &c) {
    string in, temp;
    int first, last;
    stringstream ss;

    for (int i = 0; i < c.index_1.size(); i++) {
        getline(fin, temp);
        first = temp.find('"');
        last = temp.find_last_of('"');
        in = temp.substr(first + 1, last - first - 1);
        ss << in;

        vector<double> value;
        while (getline(ss, temp, ',')) {
            value.push_back(atof(temp.c_str()));
        }
        c.rise_transition.push_back(value);

        ss.str("");
        ss.clear();
    }
}

void Parser::parse_fall_transition(fstream &fin, CellLibrary &c) {
    string in, temp;
    int first, last;
    stringstream ss;

    for (int i = 0; i < c.index_1.size(); i++) {
        getline(fin, temp);
        first = temp.find('"');
        last = temp.find_last_of('"');
        in = temp.substr(first + 1, last - first - 1);
        ss << in;

        vector<double> value;
        while (getline(ss, temp, ',')) {
            value.push_back(atof(temp.c_str()));
        }
        c.fall_transition.push_back(value);
        ss.str("");
        ss.clear();
    }
}

void Parser::parse_module(char *file, Module &m, const map<string, CellLibrary> &cells) {
    fstream fin;
    fin.open(file, ios::in);
    if (!fin) {
        cerr << "Fail to open: " << file << endl;
        exit(1);
    }

    string in, temp;
    stringstream ss;
    vector<string> tokens;

    // Parse module
    while (getline(fin, in, ';')) {
        if (in.empty()) continue;
        tokens = Parser::tokenize(in);  // Split \t \n , ( )

        if (tokens[0] == "endmodule")
            break;
        if (tokens[0] == "module") {
            m.name = tokens[1];
            continue;
        } else if (tokens[0] == "input") {
            for (int i = 1; i < tokens.size(); ++i) {
                m.input.push_back(tokens[i]);
            }
            continue;
        } else if (tokens[0] == "output") {
            for (int i = 1; i < tokens.size(); ++i) {
                m.output.push_back(tokens[i]);
            }
            continue;
        } else if (tokens[0] == "wire") {
            for (int i = 1; i < tokens.size(); ++i) {
                m.wire.push_back(tokens[i]);
            }
            continue;
        } else {    // Gates
            auto pos = cells.find(tokens[0]);
            if (pos == cells.end()) {
                cerr << "This cell " << tokens[0] << " is undefined!" << endl;
                continue;
            }
            CellLibrary cell_def = pos->second;
            Gate g;
            g.type = tokens[0];
            g.name = tokens[1];
            for (int i = 0; i < tokens.size(); ++i) {
                if (tokens[i][0] == '.') {
                    string port_name = tokens[i].substr(1);
                    string wire_name = tokens[i + 1];
                    auto port = cell_def.pins.find(port_name);
                    if (port == cell_def.pins.end()) {
                        cerr << "The port " << port_name << " is undefined!" << endl;
                    } else {
                        if (port->second.direct == "output") {
                            g.output_port = port_name;
                            g.output_wire = wire_name;
                        } else {
                            g.input_ports.push_back(port_name);
                            g.input_wires.push_back(wire_name);
                        }
                    }
                }
            }
            m.gates.push_back(g);
        }
    }
    fin.close();
}

int main(int argc, char *argv[]) {
    map<string, CellLibrary> cells;
    Module module;
    Parser parser;
    if (argc != 6) {
        cerr << "Args wrong, " << argc << " args" << endl;
        for (int i = 0; i < argc; ++i) {
            cout << argv[i] << endl;
        }
        exit(1);
    }

    char *pattern_path;
    char *lib_path;
    for (int i = 0; i < argc; ++i) {
        if (strcmp(argv[i], "-p") == 0) {
            pattern_path = argv[i + 1];
            break;
        }
        if (i == argc - 1) {
            cout << "No args -l" << endl;
            exit(1);
        }
    }
    for (int i = 0; i < argc; ++i) {
        if (strcmp(argv[i], "-l") == 0) {
            lib_path = argv[i + 1];
            break;
        }
        if (i == argc - 1) {
            cout << "No args -l" << endl;
            exit(1);
        }
    }

    parser.parse_cell_lib(lib_path, cells);
    parser.parse_module(argv[1], module, cells);
    module.topological_sort();
    module.delay(cells);
    module.output_file();
    return 0;
}