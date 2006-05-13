# Neuron max channels for Cat 15 Track 7c liberal spikes.rip
#             {neuronid: chanid}
'''
neuron2chan = {0: 21,
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
'''

# Neuron max channels for '2005-07-13, 29 - track 6 tracking, 70uV bipolar trig.tem'
#             {neuronid: chanid}
'''
neuron2chan = {0: 32,
               1: 0,
               2: 8,
               3: 52,
               4: 50,
               5: 28,
               6: 51,
               7: 27,
               8: 16}
'''

# Neuron max channels for '2005-08-02, 17 - track 5 tracking, 60uV bipolar trig.tem'
#             {neuronid: chanid}
#'''
neuron2chan = {0: 22,
               1: 22,
               2: 27,
               5: 0,
               6: 4,
               7: 47,
               8: 21,
               9: 25,
              10: 6}
#'''



# maybe should make a Polytrode class with the different designs as subclasses?

def chan2pos(polytrode='2a', chanid=None):
	if polytrode == '2a':
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
		'''
	elif polytrode == '1a':
		pass
	elif polytrode == '1b':
		pass
	elif polytrode == '1c':
		pass
	elif polytrode == '2b':
		pass
		'''
	else:
		raise ValueError, 'Unknown polytrode type', polytrode
	return SiteLoc[chanid]

print 'neuronid xcoord ycoord'
for (neuronid,chanid) in neuron2chan.iteritems():
	x,y = chan2pos(polytrode='2a', chanid=chanid) # unpack the tuple
	print neuronid, x, y

'''
Results for Cat 15 Track 7c liberal spikes.rip:

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

Results for '2005-07-13, 29 - track 6 tracking, 70uV bipolar trig.tem':

neuronid xcoord ycoord
0 28 1397
1 -28 1235
2 -28 715
3 28 1267
4 28 1072
5 28 1657
6 28 1202
7 28 1722
8 -28 195

Results for '2005-08-02, 17 - track 5 tracking, 60uV bipolar trig.tem':

neuronid xcoord ycoord
0 -28 1495
1 -28 1495
2 28 1722
5 -28 1235
6 -28 975
7 28 877
8 -28 1430
9 -28 1755
10 -28 845

'''


######################################


'''the following re's were useful in the steps for converting pascal records syntax to Python dicts of tuples syntax:
SiteLoc\[.*\]\.x =
SiteLoc\[.*\]\.y =
'''


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

NumSites := 54;
CenterX := 0;
SiteLoc[0].x := -28;
SiteLoc[0].y := 1235;
SiteLoc[1].x := -28;
SiteLoc[1].y := 1170;
SiteLoc[2].x := -28;
SiteLoc[2].y := 1105;
SiteLoc[3].x := -28;
SiteLoc[3].y := 1040;
SiteLoc[4].x := -28;
SiteLoc[4].y := 975;
SiteLoc[5].x := -28;
SiteLoc[5].y := 910;
SiteLoc[6].x := -28;
SiteLoc[6].y := 845;
SiteLoc[7].x := -28;
SiteLoc[7].y := 780;
SiteLoc[8].x := -28;
SiteLoc[8].y := 715;
SiteLoc[9].x := -28;
SiteLoc[9].y := 650;
SiteLoc[10].x := -28;
SiteLoc[10].y := 585;
SiteLoc[11].x := -28;
SiteLoc[11].y := 520;
SiteLoc[12].x := -28;
SiteLoc[12].y := 455;
SiteLoc[13].x := -28;
SiteLoc[13].y := 390;
SiteLoc[14].x := -28;
SiteLoc[14].y := 325;
SiteLoc[15].x := -28;
SiteLoc[15].y := 260;
SiteLoc[16].x := -28;
SiteLoc[16].y := 195;
SiteLoc[17].x := -28;
SiteLoc[17].y := 130;
SiteLoc[18].x := -28;
SiteLoc[18].y := 65;
SiteLoc[19].x := -28;
SiteLoc[19].y := 1300;
SiteLoc[20].x := -28;
SiteLoc[20].y := 1365;
SiteLoc[21].x := -28;
SiteLoc[21].y := 1430;
SiteLoc[22].x := -28;
SiteLoc[22].y := 1495;
SiteLoc[23].x := -28;
SiteLoc[23].y := 1560;
SiteLoc[24].x := -28;
SiteLoc[24].y := 1690;
SiteLoc[25].x := -28;
SiteLoc[25].y := 1755;
SiteLoc[26].x := -28;
SiteLoc[26].y := 1625;
SiteLoc[27].x := 28;
SiteLoc[27].y := 1722;
SiteLoc[28].x := 28;
SiteLoc[28].y := 1657;
SiteLoc[29].x := 28;
SiteLoc[29].y := 1592;
SiteLoc[30].x := 28;
SiteLoc[30].y := 1527;
SiteLoc[31].x := 28;
SiteLoc[31].y := 1462;
SiteLoc[32].x := 28;
SiteLoc[32].y := 1397;
SiteLoc[33].x := 28;
SiteLoc[33].y := 1332;
SiteLoc[34].x := 28;
SiteLoc[34].y := 32;
SiteLoc[35].x := 28;
SiteLoc[35].y := 97;
SiteLoc[36].x := 28;
SiteLoc[36].y := 162;
SiteLoc[37].x := 28;
SiteLoc[37].y := 227;
SiteLoc[38].x := 28;
SiteLoc[38].y := 292;
SiteLoc[39].x := 28;
SiteLoc[39].y := 357;
SiteLoc[40].x := 28;
SiteLoc[40].y := 422;
SiteLoc[41].x := 28;
SiteLoc[41].y := 487;
SiteLoc[42].x := 28;
SiteLoc[42].y := 552;
SiteLoc[43].x := 28;
SiteLoc[43].y := 617;
SiteLoc[44].x := 28;
SiteLoc[44].y := 682;
SiteLoc[45].x := 28;
SiteLoc[45].y := 747;
SiteLoc[46].x := 28;
SiteLoc[46].y := 812;
SiteLoc[47].x := 28;
SiteLoc[47].y := 877;
SiteLoc[48].x := 28;
SiteLoc[48].y := 942;
SiteLoc[49].x := 28;
SiteLoc[49].y := 1007;
SiteLoc[50].x := 28;
SiteLoc[50].y := 1072;
SiteLoc[51].x := 28;
SiteLoc[51].y := 1202;
SiteLoc[52].x := 28;
SiteLoc[52].y := 1267;
SiteLoc[53].x := 28;
SiteLoc[53].y := 1137;

'''
