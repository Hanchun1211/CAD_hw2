#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
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
    double cap = 0.0;
} Pin;

class Wire {
public:
    string name;
    string fan_in;
    vector<string> fan_out;
    double net_cap = 0.0;
};

class CellLibrary {
public:
    string name;
    map<string, Pin> pins;
    vector<double> index_1, index_2;
    vector<vector<double>> cell_rise;
    vector<vector<double>> cell_fall;
    vector<vector<double>> rise_transition;
    vector<vector<double>> fall_transition;
};

class Gate {
public:
    string name;
    string type;
    vector<string> input_ports;
    vector<string> input_wires;
    vector<bool> input_value;
    string output_port;
    string output_wire;
    int input_num_in_topological = 0;   // Use to compute topological sort
    bool isOutput = false;              // If this gate is an output
    vector<Gate *> next;                  // Index of following gates
    vector<Gate *> prev;                  // Index of previous gates

    bool out_value = false; // The output value of this gate
    double input_transition_delay = 0.0;    // Max of prev.gate_delay
    double output_cap = 0.0;                // Sum of next cap
    double propagation_delay = 0.0;
    double rise_fall_time = 0.0;
    double critical_gate_delay = 0.0;       // Critical path propagation_delay to this gate
    Gate *critical_prev_gate = nullptr;     // Where did critical path from

    Gate() = default;

    bool calc_output(const vector<string> &in_wires,
                     const vector<bool> &pat,
                     const map<string, CellLibrary> &cells) {
        // Compute the output value of this gate
        // Set input pattern
        input_value.clear();
        std::map<string, bool> pattern;
        if (in_wires.size() == pat.size()) {
            for (auto i = 0; i < in_wires.size(); ++i) {
                pattern[in_wires[i]] = pat[i];
            }
        } else {
            cerr << "Pattern error!" << endl;
        }

        for (auto &iw:input_wires) {
            if (pattern.find(iw) != pattern.end()) {
                input_value.push_back(pattern[iw]);
            }
        }
        if (!prev.empty()) {
            for (auto &g:prev) {
                input_value.push_back(g->out_value);
            }
        }

        // Calculate output value and input transition critical_path
        if (type == "NOR2X1") {
            if (input_value.size() != 2) {
                cerr << name << " NOR2X1 input value error!" << endl;
            } else {
                out_value = !(input_value[0] || input_value[1]);
                // Primary input critical_path is 0.0
                if (prev.size() == 2) {
                    input_transition_delay = std::max(prev[0]->rise_fall_time, prev[1]->rise_fall_time);
                } else if (prev.size() == 1) {
                    input_transition_delay = prev.front()->rise_fall_time;
                } else {
                    input_transition_delay = 0.0;
                }
            }
        } else if (type == "INVX1") {
            if (input_value.size() != 1) {
                cerr << name << " INVX1 input value error!" << endl;
            } else {
                out_value = !input_value[0];
                // Primary input critical_path is 0.0
                input_transition_delay = (prev.empty()) ? 0.0 : prev.front()->rise_fall_time;
            }
        } else if (type == "NANDX1") {
            if (input_value.size() != 2) {
                cerr << name << " NANDX1 input value error!" << endl;
            } else {
                out_value = !(input_value[0] && input_value[1]);
                // Primary input critical_path is 0.0
                if (prev.size() == 2) {
                    input_transition_delay = std::max(prev[0]->rise_fall_time, prev[1]->rise_fall_time);
                } else if (prev.size() == 1) {
                    input_transition_delay = prev.front()->rise_fall_time;
                } else {
                    input_transition_delay = 0.0;
                }
            }
        } else {
            cerr << "Undefined gate name" << endl;
        }

        return out_value;
    }
};

class Module {
public:
    string name;
    vector<string> input;
    vector<string> output;
    map<string, Wire> wires;
    vector<Gate> gates;
    vector<int> output_gates;
    vector<int> topological_order;
    double global_max_delay = 0.0;
    Gate *critical_output_gate = nullptr;
    list<string> critical_path;

    void topological_sort();

    void calc_critical_path();

    void output_file(ofstream &f);

    void calc_output_and_delay(const vector<string> &in_wires,
                               const vector<bool> &pattern,
                               map<string, CellLibrary> &cells);

private:
    static double lookup_table(double &row, double &col,
                               vector<double> &row_index, vector<double> &col_index,
                               vector<vector<double>> &table);
};

void Module::topological_sort() {
    list<int> queue;
    for (int i = 0; i < gates.size(); i++) {
        gates[i].input_num_in_topological = 0;
        gates[i].isOutput = false;

        // Find and flag output gates
        for (auto &j : output) {
            if (gates[i].output_wire == j) {
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
                    gates[i].prev.push_back(&gates[index]);
                    gates[index].next.push_back(&gates[i]);
                    if (gates[i].input_num_in_topological == gates[i].input_wires.size()) {
                        queue.push_back(i);
                    }
                    break;
                }
            }
        }
    }
}

double Module::lookup_table(double &row, double &col,
                            vector<double> &row_index,
                            vector<double> &col_index,
                            vector<vector<double>> &table) {
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

void Module::calc_critical_path() {
    global_max_delay = 0.0;
    critical_output_gate = nullptr;
    critical_path.clear();
    for (auto &idx:output_gates) {
        if (gates[idx].critical_gate_delay > global_max_delay) {
            global_max_delay = gates[idx].critical_gate_delay;
            critical_output_gate = &gates[idx];
        }
    }

    Gate *cur_gate = critical_output_gate;
    critical_path.push_front(cur_gate->output_wire);
    cur_gate = cur_gate->critical_prev_gate;
    while (cur_gate->critical_prev_gate != nullptr) {
        critical_path.push_front(cur_gate->output_wire);
        cur_gate = cur_gate->critical_prev_gate;
    }
    critical_path.push_front(cur_gate->output_wire);
    critical_path.push_front(cur_gate->input_wires.front());
}

void Module::output_file(ofstream &f) {
    f << "Longest delay = " << global_max_delay << ", the path is:" << endl;
    auto it = critical_path.begin();
    while (!critical_path.empty()) {
        f << critical_path.front();
        critical_path.pop_front();
        if (!critical_path.empty()) {
            f << " -> ";
        } else {
            f << endl << endl;
            break;
        }
    }
    for (auto &g:gates) {
        f << g.name << " " << g.out_value << " " << g.propagation_delay << " " << g.rise_fall_time << endl;
    }
    f << endl;
}

void Module::calc_output_and_delay(const vector<string> &in_wires,
                                   const vector<bool> &pattern,
                                   map<string, CellLibrary> &cells) {
    // Calculate output and critical_path for each gate by topological order
    for (auto &idx:topological_order) {
        gates[idx].calc_output(in_wires, pattern, cells);
        gates[idx].output_cap = wires[gates[idx].output_wire].net_cap;
        CellLibrary cell = cells[gates[idx].type];
        if (gates[idx].out_value) {
            gates[idx].propagation_delay = lookup_table(gates[idx].input_transition_delay, gates[idx].output_cap,
                                                        cell.index_2, cell.index_1, cell.cell_rise);
            gates[idx].rise_fall_time = lookup_table(gates[idx].input_transition_delay, gates[idx].output_cap,
                                                     cell.index_2, cell.index_1, cell.rise_transition);
        } else {
            gates[idx].propagation_delay = lookup_table(gates[idx].input_transition_delay, gates[idx].output_cap,
                                                        cell.index_2, cell.index_1, cell.cell_fall);
            gates[idx].rise_fall_time = lookup_table(gates[idx].input_transition_delay, gates[idx].output_cap,
                                                     cell.index_2, cell.index_1, cell.fall_transition);
        }

        // Calculate longest sensitize path
        string gate_type = gates[idx].type;
        if (gate_type == "INVX1") {
            if (gates[idx].prev.empty()) {   // Primary input inverter
                gates[idx].critical_prev_gate = nullptr;
                gates[idx].critical_gate_delay = gates[idx].propagation_delay;
            } else {
                gates[idx].critical_prev_gate = gates[idx].prev.front();
                gates[idx].critical_gate_delay =
                        gates[idx].propagation_delay + gates[idx].prev.front()->critical_gate_delay;
            }
        } else if (gate_type == "NOR2X1") {
            switch (gates[idx].prev.size()) {
                case 0:
                    // Primary input NOR gate
                    gates[idx].critical_prev_gate = nullptr;
                    gates[idx].critical_gate_delay = gates[idx].propagation_delay;
                    break;
                case 1:
                    // One primary input
                    if (!gates[idx].prev.front()->out_value) {
                        if (gates[idx].input_value[0] || gates[idx].input_value[1]) {
                            // Primary input is 1, don't care prev gate
                            gates[idx].critical_prev_gate = nullptr;
                            gates[idx].critical_gate_delay = gates[idx].propagation_delay;
                        } else {
                            // Both input are 0, pick prev gate path
                            gates[idx].critical_prev_gate = gates[idx].prev.front();
                            gates[idx].critical_gate_delay =
                                    gates[idx].propagation_delay + gates[idx].prev.front()->critical_gate_delay;
                        }
                    } else {
                        if (gates[idx].input_value[0] && gates[idx].input_value[1]) {
                            // Primary input is 1, don't care prev gate
                            gates[idx].critical_prev_gate = nullptr;
                            gates[idx].critical_gate_delay = gates[idx].propagation_delay;
                        } else {
                            // Primary input is 0, pick prev gate path
                            gates[idx].critical_prev_gate = gates[idx].prev.front();
                            gates[idx].critical_gate_delay =
                                    gates[idx].propagation_delay + gates[idx].prev.front()->critical_gate_delay;
                        }
                    }
                    break;
                case 2:
                    Gate *prev1, *prev2;
                    prev1 = gates[idx].prev[0];
                    prev2 = gates[idx].prev[1];
                    if (!prev1->out_value && !prev2->out_value) {
                        // Input are both 0 , pick longer path
                        Gate *longer = (prev1->critical_gate_delay > prev2->critical_gate_delay) ? prev1 : prev2;
                        gates[idx].critical_prev_gate = longer;
                        gates[idx].critical_gate_delay = gates[idx].propagation_delay + longer->critical_gate_delay;
                    } else if (prev1->out_value && prev2->out_value) {
                        // Both input are 1, pick shorter path
                        Gate *shorter = (prev1->critical_gate_delay < prev2->critical_gate_delay) ? prev1 : prev2;
                        gates[idx].critical_prev_gate = shorter;
                        gates[idx].critical_gate_delay = gates[idx].propagation_delay + shorter->critical_gate_delay;
                    } else {
                        // Only one input is 1, pick that path
                        Gate *sensitize = (prev1->out_value) ? prev1 : prev2;
                        gates[idx].critical_prev_gate = sensitize;
                        gates[idx].critical_gate_delay = gates[idx].propagation_delay + sensitize->critical_gate_delay;
                    }
                    break;
                default:
                    cerr << "NOR gate fan in error!" << endl;
            }
        } else if (gate_type == "NANDX1") {
            switch (gates[idx].prev.size()) {
                case 0:
                    // Primary input NAND gate
                    gates[idx].critical_prev_gate = nullptr;
                    gates[idx].critical_gate_delay = gates[idx].propagation_delay;
                    break;
                case 1:
                    // One primary input
                    if (gates[idx].prev.front()->out_value) {
                        if (gates[idx].input_value[0] && gates[idx].input_value[1]) {
                            // Both input are 1, pick prev gate path
                            gates[idx].critical_prev_gate = gates[idx].prev.front();
                            gates[idx].critical_gate_delay =
                                    gates[idx].propagation_delay + gates[idx].prev.front()->critical_gate_delay;
                        } else {
                            // Primary input is 0, don't care prev gate
                            gates[idx].critical_prev_gate = nullptr;
                            gates[idx].critical_gate_delay = gates[idx].propagation_delay;
                        }
                    } else {
                        if (gates[idx].input_value[0] || gates[idx].input_value[1]) {
                            // Primary input is 1, pick prev gate path
                            gates[idx].critical_prev_gate = gates[idx].prev.front();
                            gates[idx].critical_gate_delay =
                                    gates[idx].propagation_delay + gates[idx].prev.front()->critical_gate_delay;

                        } else {
                            // Primary input is 0, don't care prev gate path
                            gates[idx].critical_prev_gate = nullptr;
                            gates[idx].critical_gate_delay = gates[idx].propagation_delay;
                        }
                    }
                    break;
                case 2:
                    Gate *prev1, *prev2;
                    prev1 = gates[idx].prev[0];
                    prev2 = gates[idx].prev[1];
                    if (prev1->out_value && prev2->out_value) {
                        // Input are both 1 , pick longer path
                        Gate *longer = (prev1->critical_gate_delay > prev2->critical_gate_delay) ? prev1 : prev2;
                        gates[idx].critical_prev_gate = longer;
                        gates[idx].critical_gate_delay = gates[idx].propagation_delay + longer->critical_gate_delay;
                    } else if (!prev1->out_value && !prev2->out_value) {
                        // Both input are 0, pick shorter path
                        Gate *shorter = (prev1->critical_gate_delay < prev2->critical_gate_delay) ? prev1 : prev2;
                        gates[idx].critical_prev_gate = shorter;
                        gates[idx].critical_gate_delay = gates[idx].propagation_delay + shorter->critical_gate_delay;
                    } else {
                        // One of the input is 0, pick that path
                        Gate *sensitize = (!prev1->out_value) ? prev1 : prev2;
                        gates[idx].critical_prev_gate = sensitize;
                        gates[idx].critical_gate_delay = gates[idx].propagation_delay + sensitize->critical_gate_delay;
                    }
                    break;
                default:
                    cerr << "NAND gate fan in error!" << endl;
            }
        }
    }
}

class Parser {
public:
    vector<string> input_wires;         // Wire order of input patterns
    vector<vector<bool>> input_patterns; // Each patterns' input

    static void parse_cell_lib(char *file, map<string, CellLibrary> &cells);

    static void parse_module(char *file, Module &m, map<string, CellLibrary> &cells);

    void parse_pattern(char *pat_name);

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
    const std::regex reg(R"(\/\/[^\n\r]*|\/\*[\s\S]*\*\/|[\s\t\n]*)");    // Verilog comment format
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

    while (cnt < 4) {   // 4 timing property each cell
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
            value.push_back(strtod(temp.c_str(), nullptr));
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
            value.push_back(strtod(temp.c_str(), nullptr));
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
            value.push_back(strtod(temp.c_str(), nullptr));
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
            value.push_back(strtod(temp.c_str(), nullptr));
        }
        c.fall_transition.push_back(value);
        ss.str("");
        ss.clear();
    }
}

void Parser::parse_module(char *file, Module &m, map<string, CellLibrary> &cells) {
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
                Wire w;
                w.name = tokens[i];
                w.net_cap = 0.03;   // Primary output loading is 0.03
                m.wires.insert(pair<string, Wire>(tokens[i], w));
            }
            continue;
        } else if (tokens[0] == "wire") {
            for (int i = 1; i < tokens.size(); ++i) {
                Wire w;
                w.name = tokens[i];
                m.wires.insert(pair<string, Wire>(tokens[i], w));
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
                            m.wires[wire_name].fan_in = g.name;
                        } else {
                            g.input_ports.push_back(port_name);
                            g.input_wires.push_back(wire_name);
                            // Add wire net cap
                            m.wires[wire_name].fan_out.push_back(g.name);
                            double cap = cell_def.pins[port_name].cap;
                            m.wires[wire_name].net_cap += cap;
                        }
                    }
                }
            }
            m.gates.push_back(g);
        }
    }
    fin.close();
}

void Parser::parse_pattern(char *pat_name) {
    fstream fin;
    fin.open(pat_name, ios::in);
    if (!fin) {
        cerr << "Fail to open: " << pat_name << endl;
        exit(1);
    }

    string line;
    getline(fin, line);
    while (tokenize(line).front() != "input")
        getline(fin, line);
    // input n1, n2, n3
    // Store n1 n2 n3 in input_wires
    for (auto &w: tokenize(line)) {
        if (w == "input") continue;
        input_wires.push_back(w);
    }

    while (getline(fin, line)) {
        line = remove_comment_and_space(line);
        if (line == ".end") break;
        vector<bool> pattern;
        pattern.reserve(input_wires.size());
        for (int i = 0; i < input_wires.size(); ++i) {
            if (line[i] == '0')
                pattern.push_back(false);
            else
                pattern.push_back(true);
        }
        input_patterns.push_back(pattern);
    }
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
    // Build circuit
    parser.parse_module(argv[1], module, cells);
    module.topological_sort();
    // Input pattern
    parser.parse_pattern(pattern_path);

    // Output file
    ofstream file;
    file.open("309510145_" + module.name + ".txt", ios::out);
    if (!file) cerr << "Fail to create output file!" << endl;
    // Simulation process
    for (auto &pat:parser.input_patterns) {
        module.calc_output_and_delay(parser.input_wires, pat, cells);
        module.calc_critical_path();
        module.output_file(file);
    }

    return 0;
}