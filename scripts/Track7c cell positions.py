'''Cell positions for liberal spikes rip in Cat 15 Track7c'''
#pos = {cellid: chanid}
positions = {0: 21,
             1: 27,
             2: 20,
             3: 33,
             4: 43,
             5: 0,
             6: 20,
             8: 7,
             9: 46,
             10: 40,
             11: 49,
             12: 14,
             13: 30,
             14: 13,
             15: 9,
             16: 28,
             17: 25,
             18: 24,
             19: 39,
             20: 28,
             22: 44,
             23: 50,
             24: 41,
             25: 21,
             26: 21,
             27: 6}

#SURF ElectrodeTypes.pas definition for 2a polytrode
'''
Name := 'uMap54_2a';
Description := 'uMap54_2a, 65um spacing';
SiteSize.x := 15;
SiteSize.y := 15;
RoundSite := TRUE;
Created := FALSE;

NumPoints := 8;
Outline[0].x := -100;
Outline[0].y := 0;
Outline[1].x := -100;
Outline[1].y := 1805;
Outline[2].x := -20;
Outline[2].y := 2085;
Outline[3].x := 0;
Outline[3].y := 2185;
Outline[4].x := 20;
Outline[4].y := 2085;
Outline[5].x := 100;
Outline[5].y := 1805;
Outline[6].x := 100;
Outline[6].y := 0;
Outline[7].x := Outline[0].x;
Outline[7].y := Outline[0].y;
NumSites := 54;
CenterX := 0;
'''

# the following re was useful in the steps for converting pascal records to Python dicts of tuples: SiteLoc\[.*\]\.y =
SiteLoc = {} # init a dict
SiteLoc[0] = (-28, 1235)
SiteLoc[1] = (-28, 1170)
SiteLoc[2] = (-28, 1105)
SiteLoc[3] = (-28, 1040)
SiteLoc[4] = (-28, 975)
SiteLoc[5] = (-28, 910)
SiteLoc[6] = (-28, 845)
SiteLoc[7] = (-28, 780)
SiteLoc[8] = (-28, 715)
SiteLoc[9] = (-28, 650)
SiteLoc[10] = (-28, 585)
SiteLoc[11] = (-28, 520)
SiteLoc[12] = (-28, 455)
SiteLoc[13] = (-28, 390)
SiteLoc[14] = (-28, 325)
SiteLoc[15] = (-28, 260)
SiteLoc[16] = (-28, 195)
SiteLoc[17] = (-28, 130)
SiteLoc[18] = (-28, 65)
SiteLoc[19] = (-28, 1300)
SiteLoc[20] = (-28, 1365)
SiteLoc[21] = (-28, 1430)
SiteLoc[22] = (-28, 1495)
SiteLoc[23] = (-28, 1560)
SiteLoc[24] = (-28, 1690)
SiteLoc[25] = (-28, 1755)
SiteLoc[26] = (-28, 1625)
SiteLoc[27] = (28, 1722)
SiteLoc[28] = (28, 1657)
SiteLoc[29] = (28, 1592)
SiteLoc[30] = (28, 1527)
SiteLoc[31] = (28, 1462)
SiteLoc[32] = (28, 1397)
SiteLoc[33] = (28, 1332)
SiteLoc[34] = (28, 32)
SiteLoc[35] = (28, 97)
SiteLoc[36] = (28, 162)
SiteLoc[37] = (28, 227)
SiteLoc[38] = (28, 292)
SiteLoc[39] = (28, 357)
SiteLoc[40] = (28, 422)
SiteLoc[41] = (28, 487)
SiteLoc[42] = (28, 552)
SiteLoc[43] = (28, 617)
SiteLoc[44] = (28, 682)
SiteLoc[45] = (28, 747)
SiteLoc[46] = (28, 812)
SiteLoc[47] = (28, 877)
SiteLoc[48] = (28, 942)
SiteLoc[49] = (28, 1007)
SiteLoc[50] = (28, 1072)
SiteLoc[51] = (28, 1202)
SiteLoc[52] = (28, 1267)
SiteLoc[53] = (28, 1137)

print 'neuronid xcoord ycoord'
for (neuronid,chanid) in positions.iteritems():
	print neuronid, SiteLoc[chanid][0], SiteLoc[chanid][1]

'''
Results:

neuronid xcoord ycoord
0 -28 1430
1 28 1722
2 -28 1365
3 28 1332
4 28 617
5 -28 1235
6 -28 1365
8 -28 780
9 28 812
10 28 422
11 28 1007
12 -28 325
13 28 1527
14 -28 390
15 -28 650
16 28 1657
17 -28 1755
18 -28 1690
19 28 357
20 28 1657
22 28 682
23 28 1072
24 28 487
25 -28 1430
26 -28 1430
27 -28 845

'''