module c432 (N1, N4, N8, N11, N14, N17, N21, N24, N27, N30, N34, N37, N40,  N43, N47, N50, N53, N56, N60, N63, N66, N69, N73, N76, N79,  N86,  N89, N92, N95, N99, N102, N105, N108, N112, N115, N223, N329, N370,  N421, N430, N431, N432);
  input N1, N4, N8, N11, N14, N17, N21, N24, N27, N30, N34, N37, N40, N43, N47, N50, N53, N56, N60, N63, N66, N69, N73, N76, N79, N82, N86, N89, N92,N95, N99, N102, N105, N108, N112, N115;
  output N223, N329, N370, N421, N430, N431, N432;
  wire n188, n189, n190, n191, n192, n193, n194, n195, n196, n197, n198, n199, n200, n201, n202, n203, n204, n205, n206, n207, n208, n209, n210, n211, n212, n213, n214, n215, n216, n217, n218, n219, n220, n221, n222, n223, n224, n225, n226, n227, n228, n229, n230, n231, n232, n233, n234, n235, n236, n237, n238, n239, n240, n241, n242, n243, n244, n245, n246, n247, n248, n249, n250, n251, n252, n253, n254, n255, n256, n257, n258, n259, n260, n261, n262, n263, n264, n265, n266, n267, n268, n269, n270, n271, n272, n273, n274, n275, n276, n277, n278, n279, n280, n281, n282, n283, n284, n285, n286, n287, n288, n289, n290, n291, n292, n293, n294, n295, n296, n297, n298, n299, n300, n301, n302, n303, n304, n305, n306, n307, n308, n309, n310, n311, n312, n313, n314, n315, n316, n317, n318, n319, n320, n321, n322, n323, n324, n325, n326, n327, n328, n329, n330, n331, n332, n333, n334, n335, n336, n337, n338, n339, n340, n341, n342, n343, n344, n345, n346, n347, n348, n349, n350, n351;

  NANDX1 U189 (.A1(n188), .A2(n189), .ZN(N432));
  NANDX1 U190 (.A1(n190), .A2(n191), .ZN(n189));
  NANDX1 U191 (.A1(n192), .A2(n193), .ZN(n190));
  NANDX1 U192 (.A1(n194), .A2(n195), .ZN(n193));
  INVX1 U193 (.I(n196), .ZN(n195));
  NOR2X1 U194 (.A1(n197), .A2(n198), .ZN(n192));
  NOR2X1 U195 (.A1(n199), .A2(n200), .ZN(n197));
  NANDX1 U196 (.A1(n201), .A2(n202), .ZN(N431));
  NANDX1 U197 (.A1(n203), .A2(n204), .ZN(n202));
  INVX1 U198 (.I(n205), .ZN(n204));
  NOR2X1 U199 (.A1(n206), .A2(n207), .ZN(N421));
  NOR2X1 U200 (.A1(n208), .A2(n209), .ZN(n207));
  NANDX1 U201 (.A1(n210), .A2(n211), .ZN(n209));
  NANDX1 U202 (.A1(N8), .A2(N329), .ZN(n211));
  NANDX1 U203 (.A1(N14), .A2(N370), .ZN(n210));
  NOR2X1 U204 (.A1(n212), .A2(n213), .ZN(n206));
  NANDX1 U205 (.A1(n205), .A2(n200), .ZN(n213));
  NANDX1 U206 (.A1(n214), .A2(n215), .ZN(n200));
  NOR2X1 U207 (.A1(n216), .A2(n217), .ZN(n215));
  NOR2X1 U208 (.A1(n218), .A2(n219), .ZN(n214));
  INVX1 U209 (.I(n220), .ZN(n218));
  NANDX1 U210 (.A1(N370), .A2(N105), .ZN(n220));
  NOR2X1 U211 (.A1(n194), .A2(n199), .ZN(n205));
  NOR2X1 U212 (.A1(n221), .A2(n222), .ZN(n199));
  NANDX1 U213 (.A1(n223), .A2(n224), .ZN(n221));
  NANDX1 U214 (.A1(N86), .A2(N329), .ZN(n224));
  NANDX1 U215 (.A1(N92), .A2(N370), .ZN(n223));
  NOR2X1 U216 (.A1(n225), .A2(n226), .ZN(n194));
  NANDX1 U217 (.A1(n227), .A2(n228), .ZN(n225));
  NANDX1 U218 (.A1(N370), .A2(N79), .ZN(n228));
  NANDX1 U219 (.A1(N329), .A2(N73), .ZN(n227));
  INVX1 U220 (.I(n229), .ZN(n212));
  NOR2X1 U221 (.A1(N430), .A2(N108), .ZN(n229));
  NANDX1 U222 (.A1(n201), .A2(n203), .ZN(N430));
  NOR2X1 U223 (.A1(n198), .A2(n196), .ZN(n203));
  NOR2X1 U224 (.A1(n230), .A2(n231), .ZN(n196));
  NANDX1 U225 (.A1(n232), .A2(n233), .ZN(n230));
  NANDX1 U226 (.A1(N370), .A2(N66), .ZN(n233));
  NANDX1 U227 (.A1(N329), .A2(N60), .ZN(n232));
  NOR2X1 U228 (.A1(n234), .A2(n235), .ZN(n198));
  NANDX1 U229 (.A1(n236), .A2(n237), .ZN(n234));
  NANDX1 U230 (.A1(N47), .A2(N329), .ZN(n237));
  NANDX1 U231 (.A1(N53), .A2(N370), .ZN(n236));
  INVX1 U232 (.I(n238), .ZN(n201));
  NANDX1 U233 (.A1(n188), .A2(n191), .ZN(n238));
  NANDX1 U234 (.A1(n239), .A2(n240), .ZN(n191));
  NOR2X1 U235 (.A1(n241), .A2(n242), .ZN(n239));
  INVX1 U236 (.I(n243), .ZN(n242));
  NANDX1 U237 (.A1(N370), .A2(N40), .ZN(n243));
  NOR2X1 U238 (.A1(n244), .A2(n245), .ZN(n241));
  NANDX1 U239 (.A1(n246), .A2(n247), .ZN(n188));
  INVX1 U240 (.I(n248), .ZN(n246));
  NANDX1 U241 (.A1(n249), .A2(n250), .ZN(n248));
  NANDX1 U242 (.A1(N21), .A2(N329), .ZN(n250));
  NANDX1 U243 (.A1(N27), .A2(N370), .ZN(n249));
  NANDX1 U244 (.A1(n251), .A2(n252), .ZN(N370));
  NOR2X1 U245 (.A1(n253), .A2(n254), .ZN(n252));
  NANDX1 U246 (.A1(n255), .A2(n256), .ZN(n254));
  NANDX1 U247 (.A1(n257), .A2(n247), .ZN(n256));
  INVX1 U248 (.I(n258), .ZN(n247));
  NOR2X1 U249 (.A1(N27), .A2(n259), .ZN(n257));
  NOR2X1 U250 (.A1(n260), .A2(n244), .ZN(n259));
  NANDX1 U251 (.A1(n261), .A2(n262), .ZN(n255));
  INVX1 U252 (.I(n222), .ZN(n262));
  NOR2X1 U253 (.A1(N92), .A2(n263), .ZN(n261));
  NOR2X1 U254 (.A1(n264), .A2(n244), .ZN(n263));
  NANDX1 U255 (.A1(n265), .A2(n266), .ZN(n253));
  NANDX1 U256 (.A1(n267), .A2(n268), .ZN(n266));
  INVX1 U257 (.I(n208), .ZN(n268));
  NOR2X1 U258 (.A1(N14), .A2(n269), .ZN(n267));
  NOR2X1 U259 (.A1(n270), .A2(n244), .ZN(n269));
  NOR2X1 U260 (.A1(n271), .A2(n272), .ZN(n265));
  NOR2X1 U261 (.A1(n235), .A2(n273), .ZN(n272));
  NANDX1 U262 (.A1(n274), .A2(n275), .ZN(n273));
  INVX1 U263 (.I(N53), .ZN(n275));
  NANDX1 U264 (.A1(N329), .A2(n276), .ZN(n274));
  NOR2X1 U265 (.A1(n277), .A2(n278), .ZN(n271));
  NANDX1 U266 (.A1(n279), .A2(n280), .ZN(n278));
  INVX1 U267 (.I(N40), .ZN(n280));
  NANDX1 U268 (.A1(N329), .A2(n281), .ZN(n279));
  NOR2X1 U269 (.A1(n282), .A2(n283), .ZN(n251));
  NANDX1 U270 (.A1(n284), .A2(n285), .ZN(n283));
  NANDX1 U271 (.A1(n286), .A2(n287), .ZN(n285));
  NOR2X1 U272 (.A1(N105), .A2(n216), .ZN(n287));
  NOR2X1 U273 (.A1(n217), .A2(n219), .ZN(n286));
  INVX1 U274 (.I(N95), .ZN(n219));
  INVX1 U275 (.I(n288), .ZN(n217));
  NANDX1 U276 (.A1(N99), .A2(N329), .ZN(n288));
  NANDX1 U277 (.A1(n289), .A2(n290), .ZN(n284));
  NOR2X1 U278 (.A1(N115), .A2(n291), .ZN(n290));
  NOR2X1 U279 (.A1(n292), .A2(n293), .ZN(n289));
  NOR2X1 U280 (.A1(n294), .A2(n244), .ZN(n292));
  INVX1 U281 (.I(n295), .ZN(n294));
  NANDX1 U282 (.A1(n296), .A2(n297), .ZN(n282));
  NANDX1 U283 (.A1(n298), .A2(n299), .ZN(n297));
  INVX1 U284 (.I(n231), .ZN(n299));
  NOR2X1 U285 (.A1(N66), .A2(n300), .ZN(n298));
  NOR2X1 U286 (.A1(n301), .A2(n244), .ZN(n300));
  NANDX1 U287 (.A1(n302), .A2(n303), .ZN(n296));
  INVX1 U288 (.I(n226), .ZN(n303));
  NOR2X1 U289 (.A1(N79), .A2(n304), .ZN(n302));
  NOR2X1 U290 (.A1(n305), .A2(n244), .ZN(n304));
  INVX1 U291 (.I(N329), .ZN(n244));
  NANDX1 U292 (.A1(n306), .A2(n307), .ZN(N329));
  NOR2X1 U293 (.A1(n308), .A2(n309), .ZN(n307));
  INVX1 U294 (.I(n310), .ZN(n309));
  NOR2X1 U295 (.A1(n260), .A2(n270), .ZN(n310));
  NOR2X1 U296 (.A1(n208), .A2(N8), .ZN(n270));
  NANDX1 U297 (.A1(N4), .A2(n311), .ZN(n208));
  NANDX1 U298 (.A1(N1), .A2(N223), .ZN(n311));
  NOR2X1 U299 (.A1(n258), .A2(N21), .ZN(n260));
  NANDX1 U300 (.A1(N17), .A2(n312), .ZN(n258));
  NANDX1 U301 (.A1(N11), .A2(N223), .ZN(n312));
  NANDX1 U302 (.A1(n313), .A2(n295), .ZN(n308));
  NANDX1 U303 (.A1(n314), .A2(N108), .ZN(n295));
  NOR2X1 U304 (.A1(N112), .A2(n291), .ZN(n314));
  INVX1 U305 (.I(n315), .ZN(n291));
  NANDX1 U306 (.A1(N102), .A2(N223), .ZN(n315));
  NOR2X1 U307 (.A1(n264), .A2(n305), .ZN(n313));
  NOR2X1 U308 (.A1(n226), .A2(N73), .ZN(n305));
  NANDX1 U309 (.A1(N69), .A2(n316), .ZN(n226));
  NANDX1 U310 (.A1(N63), .A2(N223), .ZN(n316));
  NOR2X1 U311 (.A1(n222), .A2(N86), .ZN(n264));
  NANDX1 U312 (.A1(N82), .A2(n317), .ZN(n222));
  NANDX1 U313 (.A1(N76), .A2(N223), .ZN(n317));
  NOR2X1 U314 (.A1(n318), .A2(n319), .ZN(n306));
  NANDX1 U315 (.A1(n320), .A2(n276), .ZN(n319));
  INVX1 U316 (.I(n321), .ZN(n276));
  NOR2X1 U317 (.A1(n235), .A2(N47), .ZN(n321));
  NANDX1 U318 (.A1(N43), .A2(n322), .ZN(n235));
  NANDX1 U319 (.A1(N37), .A2(N223), .ZN(n322));
  NANDX1 U320 (.A1(n323), .A2(N95), .ZN(n320));
  NOR2X1 U321 (.A1(N99), .A2(n216), .ZN(n323));
  NOR2X1 U322 (.A1(n324), .A2(n325), .ZN(n216));
  INVX1 U323 (.I(N223), .ZN(n325));
  NANDX1 U324 (.A1(n326), .A2(n281), .ZN(n318));
  NANDX1 U325 (.A1(n240), .A2(n245), .ZN(n281));
  INVX1 U326 (.I(N34), .ZN(n245));
  INVX1 U327 (.I(n277), .ZN(n240));
  NANDX1 U328 (.A1(N30), .A2(n327), .ZN(n277));
  NANDX1 U329 (.A1(N24), .A2(N223), .ZN(n327));
  INVX1 U330 (.I(n301), .ZN(n326));
  NOR2X1 U331 (.A1(n231), .A2(N60), .ZN(n301));
  NANDX1 U332 (.A1(N56), .A2(n328), .ZN(n231));
  NANDX1 U333 (.A1(N50), .A2(N223), .ZN(n328));
  NANDX1 U334 (.A1(n329), .A2(n330), .ZN(N223));
  NOR2X1 U335 (.A1(n331), .A2(n332), .ZN(n330));
  NANDX1 U336 (.A1(n333), .A2(n334), .ZN(n332));
  NANDX1 U337 (.A1(N30), .A2(n335), .ZN(n334));
  INVX1 U338 (.I(N24), .ZN(n335));
  NANDX1 U339 (.A1(N43), .A2(n336), .ZN(n333));
  INVX1 U340 (.I(N37), .ZN(n336));
  NANDX1 U341 (.A1(n337), .A2(n338), .ZN(n331));
  NANDX1 U342 (.A1(N17), .A2(n339), .ZN(n338));
  INVX1 U343 (.I(N11), .ZN(n339));
  NOR2X1 U344 (.A1(n340), .A2(n341), .ZN(n337));
  NOR2X1 U345 (.A1(N102), .A2(n293), .ZN(n341));
  INVX1 U346 (.I(N108), .ZN(n293));
  NOR2X1 U347 (.A1(N1), .A2(n342), .ZN(n340));
  INVX1 U348 (.I(N4), .ZN(n342));
  NOR2X1 U349 (.A1(n343), .A2(n344), .ZN(n329));
  NANDX1 U350 (.A1(n345), .A2(n346), .ZN(n344));
  NANDX1 U351 (.A1(N82), .A2(n347), .ZN(n346));
  INVX1 U352 (.I(N76), .ZN(n347));
  NANDX1 U353 (.A1(N95), .A2(n324), .ZN(n345));
  INVX1 U354 (.I(N89), .ZN(n324));
  NANDX1 U355 (.A1(n348), .A2(n349), .ZN(n343));
  NANDX1 U356 (.A1(N56), .A2(n350), .ZN(n349));
  INVX1 U357 (.I(N50), .ZN(n350));
  NANDX1 U358 (.A1(N69), .A2(n351), .ZN(n348));
  INVX1 U359 (.I(N63), .ZN(n351));
endmodule