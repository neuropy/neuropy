# Neuron max channels for 2a polytrode Cat 15 Track 7c liberal spikes.rip
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


ypos2neuron = [(325, 12), (357, 19), (390, 14), (422, 10), (487, 24), (617, 4), (650, 15), (682, 22), (780, 8), (812, 9), (845, 27), (1007, 11), (1072, 23), (1235, 5), (1332, 3), (1365, 2), (1365, 6), (1430, 0), (1430, 25), (1430, 26), (1527, 13), (1657, 16), (1657, 20), (1690, 18), (1722, 1), (1755, 17)]

'''

# Neuron max channels for 2a polytrode Cat 15 Track 6 '2005-07-13, 29 - track 6 tracking, 70uV bipolar trig.tem'
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

# Neuron max channels for 2a polytrode Cat 15 Track 5 '2005-08-02, 17 - track 5 tracking, 60uV bipolar trig.tem'
#             {neuronid: chanid}
'''
neuron2chan = {0: 22,
               1: 22,
               2: 27,
               5: 0,
               6: 4,
               7: 47,
               8: 21,
               9: 25,
              10: 6}
'''

# Neuron max channels for 1b polytrode Cat 11 Track 1 '2004-05-07, 11 - +-49uV trig, whole file random sample.bak.tem; 2006-05-14 rip, some templates deleted by MAS.tem'
#             {neuronid: chanid}
'''
neuron2chan = {0: 46,
               1: 11,
               2: 14,
               3: 42,
               4: 53,
               5: 17,
               6: 51,
               7: 31,
               8: 15,
               9: 24,
              10: 35,
              11: 19,
              12: 22,
              13: 22,
              14: 2,
              15: 29,
              16: 13,
              17: 36,
              18: 27,
              19: 41,
              20: 5,
              21: 31,
              22: 27,
              23: 45,
              24: 24,
              25: 19,
              26: 50,
              27: 49,
              28: 0,
              29: 10,
              30: 50,
              31: 25,
              32: 4,
              33: 25,
              34: 39,
              35: 0,
              36: 10,
              37: 47,
              38: 17}
'''


# Neuron max channels for 2b polytrode Cat 13 Track 8 '2004-04-28 - 07 - +-55uV trig, 0-10 mins random sample.tem; 2006-05-15 rip, some templates deleted by MAS.tem'
#             {neuronid: chanid}
#'''
neuron2chan = {0: 29,
               1: 28,
               2: 23,
               3: 22,
               4: 3,
               5: 0,
               6: 49,
               7: 5,
               8: 53,
               9: 52,
              10: 24,
              11: 13,
              12: 24,
              13: 12,
              14: 13,
              15: 12,
              16: 6,
              17: 23,
              18: 25,
              19: 53,
              20: 30,
              21: 14,
              22: 15,
              23: 17,
              24: 11,
              25: 25,
              26: 41,
              27: 10,
              28: 27,
              29: 11}
#'''


"""Here's a BaseLayout class with the different Polytrode designs as subclasses"""

class BaseLayout(object):
    pass


class Polytrode_1a(BaseLayout):
    def __init__(self):
        self.layout = '1a'
        self.description = 'uMap54_1a, 65um spacing, 3 column hexagonal'
        SiteLoc = {}
        SiteLoc[0] = (-56, 1170)
        SiteLoc[1] = (-56, 1105)
        SiteLoc[2] = (-56, 1040)
        SiteLoc[3] = (-56, 975)
        SiteLoc[4] = (-56, 910)
        SiteLoc[5] = (-56, 845)
        SiteLoc[6] = (-56, 585)
        SiteLoc[7] = (-56, 455)
        SiteLoc[8] = (-56, 325)
        SiteLoc[9] = (-56, 195)
        SiteLoc[10] = (-56, 65)
        SiteLoc[11] = (-56, 130)
        SiteLoc[12] = (-56, 260)
        SiteLoc[13] = (-56, 390)
        SiteLoc[14] = (-56, 520)
        SiteLoc[15] = (-56, 650)
        SiteLoc[16] = (-56, 715)
        SiteLoc[17] = (-56, 780)
        SiteLoc[18] = (0, 1072)
        SiteLoc[19] = (0, 942)
        SiteLoc[20] = (0, 812)
        SiteLoc[21] = (0, 682)
        SiteLoc[22] = (0, 552)
        SiteLoc[23] = (0, 422)
        SiteLoc[24] = (0, 162)
        SiteLoc[25] = (0, 97)
        SiteLoc[26] = (0, 292)
        SiteLoc[27] = (0, 227)
        SiteLoc[28] = (0, 357)
        SiteLoc[29] = (0, 487)
        SiteLoc[30] = (0, 617)
        SiteLoc[31] = (0, 747)
        SiteLoc[32] = (0, 877)
        SiteLoc[33] = (0, 1007)
        SiteLoc[34] = (0, 1137)
        SiteLoc[35] = (0, 1202)
        SiteLoc[36] = (56, 780)
        SiteLoc[37] = (56, 650)
        SiteLoc[38] = (56, 520)
        SiteLoc[39] = (56, 390)
        SiteLoc[40] = (56, 260)
        SiteLoc[41] = (56, 130)
        SiteLoc[42] = (56, 65)
        SiteLoc[43] = (56, 195)
        SiteLoc[44] = (56, 325)
        SiteLoc[45] = (56, 455)
        SiteLoc[46] = (56, 585)
        SiteLoc[47] = (56, 715)
        SiteLoc[48] = (56, 845)
        SiteLoc[49] = (56, 910)
        SiteLoc[50] = (56, 975)
        SiteLoc[51] = (56, 1105)
        SiteLoc[52] = (56, 1170)
        SiteLoc[53] = (56, 1040)
        self.SiteLoc = SiteLoc


class Polytrode_1b(BaseLayout):
    def __init__(self):
        self.layout = '1b'
        self.description = 'uMap54_1b, 50um horzontal/46um vertical spacing, 3 column collinear'
        SiteLoc = {}
        SiteLoc[0] = (-43, 900)
        SiteLoc[1] = (-43, 850)
        SiteLoc[2] = (-43, 800)
        SiteLoc[3] = (-43, 750)
        SiteLoc[4] = (-43, 700)
        SiteLoc[5] = (-43, 650)
        SiteLoc[6] = (-43, 600)
        SiteLoc[7] = (-43, 550)
        SiteLoc[8] = (-43, 500)
        SiteLoc[9] = (-43, 450)
        SiteLoc[10] = (-43, 400)
        SiteLoc[11] = (-43, 350)
        SiteLoc[12] = (-43, 300)
        SiteLoc[13] = (-43, 250)
        SiteLoc[14] = (-43, 200)
        SiteLoc[15] = (-43, 150)
        SiteLoc[16] = (-43, 50)
        SiteLoc[17] = (-43, 100)
        SiteLoc[18] = (0, 900)
        SiteLoc[19] = (0, 800)
        SiteLoc[20] = (0, 700)
        SiteLoc[21] = (0, 600)
        SiteLoc[22] = (0, 500)
        SiteLoc[23] = (0, 400)
        SiteLoc[24] = (0, 200)
        SiteLoc[25] = (0, 100)
        SiteLoc[26] = (0, 300)
        SiteLoc[27] = (0, 50)
        SiteLoc[28] = (0, 150)
        SiteLoc[29] = (0, 250)
        SiteLoc[30] = (0, 350)
        SiteLoc[31] = (0, 450)
        SiteLoc[32] = (0, 550)
        SiteLoc[33] = (0, 650)
        SiteLoc[34] = (0, 750)
        SiteLoc[35] = (0, 850)
        SiteLoc[36] = (43, 200)
        SiteLoc[37] = (43, 100)
        SiteLoc[38] = (43, 50)
        SiteLoc[39] = (43, 150)
        SiteLoc[40] = (43, 250)
        SiteLoc[41] = (43, 300)
        SiteLoc[42] = (43, 350)
        SiteLoc[43] = (43, 400)
        SiteLoc[44] = (43, 450)
        SiteLoc[45] = (43, 500)
        SiteLoc[46] = (43, 550)
        SiteLoc[47] = (43, 600)
        SiteLoc[48] = (43, 650)
        SiteLoc[49] = (43, 700)
        SiteLoc[50] = (43, 750)
        SiteLoc[51] = (43, 850)
        SiteLoc[52] = (43, 900)
        SiteLoc[53] = (43, 800)
        self.SiteLoc = SiteLoc


class Polytrode_1c(BaseLayout):
    def __init__(self):
        self.layout = '1c'
        self.description = 'uMap54_1c, 75um spacing, 3 column, hexagonal'
        SiteLoc = {}
        SiteLoc[0] = (-65, 1251)
        SiteLoc[1] = (-65, 1101)
        SiteLoc[2] = (-65, 951)
        SiteLoc[3] = (-65, 801)
        SiteLoc[4] = (-65, 651)
        SiteLoc[5] = (-65, 501)
        SiteLoc[6] = (-65, 351)
        SiteLoc[7] = (-65, 201)
        SiteLoc[8] = (-65, 51)
        SiteLoc[9] = (-65, 126)
        SiteLoc[10] = (-65, 276)
        SiteLoc[11] = (-65, 426)
        SiteLoc[12] = (-65, 576)
        SiteLoc[13] = (-65, 726)
        SiteLoc[14] = (-65, 876)
        SiteLoc[15] = (-65, 1026)
        SiteLoc[16] = (-65, 1176)
        SiteLoc[17] = (-65, 1326)
        SiteLoc[18] = (0, 1364)
        SiteLoc[19] = (0, 1214)
        SiteLoc[20] = (0, 1064)
        SiteLoc[21] = (0, 914)
        SiteLoc[22] = (0, 764)
        SiteLoc[23] = (0, 614)
        SiteLoc[24] = (0, 314)
        SiteLoc[25] = (0, 164)
        SiteLoc[26] = (0, 464)
        SiteLoc[27] = (0, 89)
        SiteLoc[28] = (0, 239)
        SiteLoc[29] = (0, 389)
        SiteLoc[30] = (0, 539)
        SiteLoc[31] = (0, 689)
        SiteLoc[32] = (0, 839)
        SiteLoc[33] = (0, 989)
        SiteLoc[34] = (0, 1139)
        SiteLoc[35] = (0, 1289)
        SiteLoc[36] = (65, 1326)
        SiteLoc[37] = (65, 1251)
        SiteLoc[38] = (65, 1176)
        SiteLoc[39] = (65, 1026)
        SiteLoc[40] = (65, 876)
        SiteLoc[41] = (65, 726)
        SiteLoc[42] = (65, 576)
        SiteLoc[43] = (65, 426)
        SiteLoc[44] = (65, 276)
        SiteLoc[45] = (65, 126)
        SiteLoc[46] = (65, 51)
        SiteLoc[47] = (65, 201)
        SiteLoc[48] = (65, 351)
        SiteLoc[49] = (65, 501)
        SiteLoc[50] = (65, 651)
        SiteLoc[51] = (65, 951)
        SiteLoc[52] = (65, 1101)
        SiteLoc[53] = (65, 801)
        self.SiteLoc = SiteLoc


class Polytrode_2a(BaseLayout):
    def __init__(self):
        self.layout = '2a'
        self.description = 'uMap54_2a, 65um spacing, 2 column, staggered'
        SiteLoc = {}
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
        self.SiteLoc = SiteLoc


class Polytrode_2b(BaseLayout):
    def __init__(self):
        self.layout = '2b'
        self.description = 'uMap54_2b, 50um spacing, 2 column, staggered'
        SiteLoc = {}
        SiteLoc[0] = (-25, 1275)
        SiteLoc[1] = (-25, 1175)
        SiteLoc[2] = (-25, 1075)
        SiteLoc[3] = (-25, 975)
        SiteLoc[4] = (-25, 875)
        SiteLoc[5] = (-25, 775)
        SiteLoc[6] = (-25, 725)
        SiteLoc[7] = (-25, 675)
        SiteLoc[8] = (-25, 625)
        SiteLoc[9] = (-25, 575)
        SiteLoc[10] = (-25, 525)
        SiteLoc[11] = (-25, 475)
        SiteLoc[12] = (-25, 425)
        SiteLoc[13] = (-25, 375)
        SiteLoc[14] = (-25, 325)
        SiteLoc[15] = (-25, 275)
        SiteLoc[16] = (-25, 225)
        SiteLoc[17] = (-25, 175)
        SiteLoc[18] = (-25, 125)
        SiteLoc[19] = (-25, 75)
        SiteLoc[20] = (-25, 25)
        SiteLoc[21] = (-25, 825)
        SiteLoc[22] = (-25, 925)
        SiteLoc[23] = (-25, 1025)
        SiteLoc[24] = (-25, 1225)
        SiteLoc[25] = (-25, 1325)
        SiteLoc[26] = (-25, 1125)
        SiteLoc[27] = (25, 1300)
        SiteLoc[28] = (25, 1200)
        SiteLoc[29] = (25, 1100)
        SiteLoc[30] = (25, 1000)
        SiteLoc[31] = (25, 900)
        SiteLoc[32] = (25, 0)
        SiteLoc[33] = (25, 50)
        SiteLoc[34] = (25, 100)
        SiteLoc[35] = (25, 150)
        SiteLoc[36] = (25, 200)
        SiteLoc[37] = (25, 250)
        SiteLoc[38] = (25, 300)
        SiteLoc[39] = (25, 350)
        SiteLoc[40] = (25, 400)
        SiteLoc[41] = (25, 450)
        SiteLoc[42] = (25, 500)
        SiteLoc[43] = (25, 550)
        SiteLoc[44] = (25, 600)
        SiteLoc[45] = (25, 650)
        SiteLoc[46] = (25, 700)
        SiteLoc[47] = (25, 750)
        SiteLoc[48] = (25, 800)
        SiteLoc[49] = (25, 850)
        SiteLoc[50] = (25, 950)
        SiteLoc[51] = (25, 1150)
        SiteLoc[52] = (25, 1250)
        SiteLoc[53] = (25, 1050)
        self.SiteLoc = SiteLoc


#################################################################


def chan2pos(polytrode=None, chanid=None):
    SiteLoc = {} # init a dict
    if polytrode == '1a':
        """Need to fill it in"""
        pass
    elif polytrode == '1b':
        SiteLoc[0] = (-43, 900)
        SiteLoc[1] = (-43, 850)
        SiteLoc[2] = (-43, 800)
        SiteLoc[3] = (-43, 750)
        SiteLoc[4] = (-43, 700)
        SiteLoc[5] = (-43, 650)
        SiteLoc[6] = (-43, 600)
        SiteLoc[7] = (-43, 550)
        SiteLoc[8] = (-43, 500)
        SiteLoc[9] = (-43, 450)
        SiteLoc[10] = (-43, 400)
        SiteLoc[11] = (-43, 350)
        SiteLoc[12] = (-43, 300)
        SiteLoc[13] = (-43, 250)
        SiteLoc[14] = (-43, 200)
        SiteLoc[15] = (-43, 150)
        SiteLoc[16] = (-43, 50)
        SiteLoc[17] = (-43, 100)
        SiteLoc[18] = (0, 900)
        SiteLoc[19] = (0, 800)
        SiteLoc[20] = (0, 700)
        SiteLoc[21] = (0, 600)
        SiteLoc[22] = (0, 500)
        SiteLoc[23] = (0, 400)
        SiteLoc[24] = (0, 200)
        SiteLoc[25] = (0, 100)
        SiteLoc[26] = (0, 300)
        SiteLoc[27] = (0, 50)
        SiteLoc[28] = (0, 150)
        SiteLoc[29] = (0, 250)
        SiteLoc[30] = (0, 350)
        SiteLoc[31] = (0, 450)
        SiteLoc[32] = (0, 550)
        SiteLoc[33] = (0, 650)
        SiteLoc[34] = (0, 750)
        SiteLoc[35] = (0, 850)
        SiteLoc[36] = (43, 200)
        SiteLoc[37] = (43, 100)
        SiteLoc[38] = (43, 50)
        SiteLoc[39] = (43, 150)
        SiteLoc[40] = (43, 250)
        SiteLoc[41] = (43, 300)
        SiteLoc[42] = (43, 350)
        SiteLoc[43] = (43, 400)
        SiteLoc[44] = (43, 450)
        SiteLoc[45] = (43, 500)
        SiteLoc[46] = (43, 550)
        SiteLoc[47] = (43, 600)
        SiteLoc[48] = (43, 650)
        SiteLoc[49] = (43, 700)
        SiteLoc[50] = (43, 750)
        SiteLoc[51] = (43, 850)
        SiteLoc[52] = (43, 900)
        SiteLoc[53] = (43, 800)
    elif polytrode == '1c':
        """Need to fill it in"""
        pass
    elif polytrode == '2a':
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
    elif polytrode == '2b':
        SiteLoc[0] = (-25, 1275)
        SiteLoc[1] = (-25, 1175)
        SiteLoc[2] = (-25, 1075)
        SiteLoc[3] = (-25, 975)
        SiteLoc[4] = (-25, 875)
        SiteLoc[5] = (-25, 775)
        SiteLoc[6] = (-25, 725)
        SiteLoc[7] = (-25, 675)
        SiteLoc[8] = (-25, 625)
        SiteLoc[9] = (-25, 575)
        SiteLoc[10] = (-25, 525)
        SiteLoc[11] = (-25, 475)
        SiteLoc[12] = (-25, 425)
        SiteLoc[13] = (-25, 375)
        SiteLoc[14] = (-25, 325)
        SiteLoc[15] = (-25, 275)
        SiteLoc[16] = (-25, 225)
        SiteLoc[17] = (-25, 175)
        SiteLoc[18] = (-25, 125)
        SiteLoc[19] = (-25, 75)
        SiteLoc[20] = (-25, 25)
        SiteLoc[21] = (-25, 825)
        SiteLoc[22] = (-25, 925)
        SiteLoc[23] = (-25, 1025)
        SiteLoc[24] = (-25, 1225)
        SiteLoc[25] = (-25, 1325)
        SiteLoc[26] = (-25, 1125)
        SiteLoc[27] = (25, 1300)
        SiteLoc[28] = (25, 1200)
        SiteLoc[29] = (25, 1100)
        SiteLoc[30] = (25, 1000)
        SiteLoc[31] = (25, 900)
        SiteLoc[32] = (25, 0)
        SiteLoc[33] = (25, 50)
        SiteLoc[34] = (25, 100)
        SiteLoc[35] = (25, 150)
        SiteLoc[36] = (25, 200)
        SiteLoc[37] = (25, 250)
        SiteLoc[38] = (25, 300)
        SiteLoc[39] = (25, 350)
        SiteLoc[40] = (25, 400)
        SiteLoc[41] = (25, 450)
        SiteLoc[42] = (25, 500)
        SiteLoc[43] = (25, 550)
        SiteLoc[44] = (25, 600)
        SiteLoc[45] = (25, 650)
        SiteLoc[46] = (25, 700)
        SiteLoc[47] = (25, 750)
        SiteLoc[48] = (25, 800)
        SiteLoc[49] = (25, 850)
        SiteLoc[50] = (25, 950)
        SiteLoc[51] = (25, 1150)
        SiteLoc[52] = (25, 1250)
        SiteLoc[53] = (25, 1050)
    else:
        raise ValueError, 'Unknown polytrode type', polytrode
    try:
        return SiteLoc[chanid] # return the (x, y) tuple
    except TypeError: # no chanid was specified, return the whole SiteLoc dict
        return SiteLoc


print 'neuronid xcoord ycoord'
for neuronid, chanid in neuron2chan.iteritems():
    x, y = chan2pos(polytrode='2b', chanid=chanid) # unpack the tuple
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

Results for '2004-05-07, 11 - +-49uV trig, whole file random sample.bak.tem; 2006-05-14 rip, some templates deleted by MAS.tem':

neuronid xcoord ycoord
0 43 550
1 -43 350
2 -43 200
3 43 350
4 43 800
5 -43 100
6 43 850
7 0 450
8 -43 150
9 0 200
10 0 850
11 0 800
12 0 500
13 0 500
14 -43 800
15 0 250
16 -43 250
17 43 200
18 0 50
19 43 300
20 -43 650
21 0 450
22 0 50
23 43 500
24 0 200
25 0 800
26 43 750
27 43 700
28 -43 900
29 -43 400
30 43 750
31 0 100
32 -43 700
33 0 100
34 43 150
35 -43 900
36 -43 400
37 43 600
38 -43 100

Results for 2b polytrode Cat 13 Track 8 '2004-04-28 - 07 - +-55uV trig, 0-10 mins random sample.tem; 2006-05-15 rip, some templates deleted by MAS.tem':

neuronid xcoord ycoord
0 25 1100
1 25 1200
2 -25 1025
3 -25 925
4 -25 975
5 -25 1275
6 25 850
7 -25 775
8 25 1050
9 25 1250
10 -25 1225
11 -25 375
12 -25 1225
13 -25 425
14 -25 375
15 -25 425
16 -25 725
17 -25 1025
18 -25 1325
19 25 1050
20 25 1000
21 -25 325
22 -25 275
23 -25 175
24 -25 475
25 -25 1325
26 25 450
27 -25 525
28 25 1300
29 -25 475

'''


######################################


"""the following re's were useful in the steps for converting pascal records syntax to Python dicts of tuples syntax:

1)  '.x :=' replaced with ' = ('
2)  ';\n.*SiteLoc\[.*\]\.y :=' replaced with ','
3)  ';' replaced with ')'

SiteLoc\[.*\]\.x :=
SiteLoc\[.*\]\.y :=
"""
#SURF ElectrodeTypes.pas definition for 1b polytrode
'''
Name := 'uMap54_1b';
Description := 'uMap54_1b, 50um spacing';
SiteSize.x := 15;
SiteSize.y := 15;
RoundSite := TRUE;
Created := FALSE;

NumPoints := 8;
Outline[0].x := -100;
Outline[0].y := 0;
Outline[1].x := -100;
Outline[1].y := 900;
Outline[2].x := -20;
Outline[2].y := 1200;
Outline[3].x := 0;
Outline[3].y := 1250;
Outline[4].x := 20;
Outline[4].y := 1200;
Outline[5].x := 100;
Outline[5].y := 900;
Outline[6].x := 100;
Outline[6].y := 0;
Outline[7].x := Outline[0].x;
Outline[7].y := Outline[0].y;

NumSites := 54;
CenterX := 0;
SiteLoc[0].x := -43;
SiteLoc[0].y := 900;
SiteLoc[1].x := -43;
SiteLoc[1].y := 850;
SiteLoc[2].x := -43;
SiteLoc[2].y := 800;
SiteLoc[3].x := -43;
SiteLoc[3].y := 750;
SiteLoc[4].x := -43;
SiteLoc[4].y := 700;
SiteLoc[5].x := -43;
SiteLoc[5].y := 650;
SiteLoc[6].x := -43;
SiteLoc[6].y := 600;
SiteLoc[7].x := -43;
SiteLoc[7].y := 550;
SiteLoc[8].x := -43;
SiteLoc[8].y := 500;
SiteLoc[9].x := -43;
SiteLoc[9].y := 450;
SiteLoc[10].x := -43;
SiteLoc[10].y := 400;
SiteLoc[11].x := -43;
SiteLoc[11].y := 350;
SiteLoc[12].x := -43;
SiteLoc[12].y := 300;
SiteLoc[13].x := -43;
SiteLoc[13].y := 250;
SiteLoc[14].x := -43;
SiteLoc[14].y := 200;
SiteLoc[15].x := -43;
SiteLoc[15].y := 150;
SiteLoc[16].x := -43;
SiteLoc[16].y := 50;
SiteLoc[17].x := -43;
SiteLoc[17].y := 100;
SiteLoc[18].x := 0;
SiteLoc[18].y := 900;
SiteLoc[19].x := 0;
SiteLoc[19].y := 800;
SiteLoc[20].x := 0;
SiteLoc[20].y := 700;
SiteLoc[21].x := 0;
SiteLoc[21].y := 600;
SiteLoc[22].x := 0;
SiteLoc[22].y := 500;
SiteLoc[23].x := 0;
SiteLoc[23].y := 400;
SiteLoc[24].x := 0;
SiteLoc[24].y := 200;
SiteLoc[25].x := 0;
SiteLoc[25].y := 100;
SiteLoc[26].x := 0;
SiteLoc[26].y := 300;
SiteLoc[27].x := 0;
SiteLoc[27].y := 50;
SiteLoc[28].x := 0;
SiteLoc[28].y := 150;
SiteLoc[29].x := 0;
SiteLoc[29].y := 250;
SiteLoc[30].x := 0;
SiteLoc[30].y := 350;
SiteLoc[31].x := 0;
SiteLoc[31].y := 450;
SiteLoc[32].x := 0;
SiteLoc[32].y := 550;
SiteLoc[33].x := 0;
SiteLoc[33].y := 650;
SiteLoc[34].x := 0;
SiteLoc[34].y := 750;
SiteLoc[35].x := 0;
SiteLoc[35].y := 850;
SiteLoc[36].x := 43;
SiteLoc[36].y := 200;
SiteLoc[37].x := 43;
SiteLoc[37].y := 100;
SiteLoc[38].x := 43;
SiteLoc[38].y := 50;
SiteLoc[39].x := 43;
SiteLoc[39].y := 150;
SiteLoc[40].x := 43;
SiteLoc[40].y := 250;
SiteLoc[41].x := 43;
SiteLoc[41].y := 300;
SiteLoc[42].x := 43;
SiteLoc[42].y := 350;
SiteLoc[43].x := 43;
SiteLoc[43].y := 400;
SiteLoc[44].x := 43;
SiteLoc[44].y := 450;
SiteLoc[45].x := 43;
SiteLoc[45].y := 500;
SiteLoc[46].x := 43;
SiteLoc[46].y := 550;
SiteLoc[47].x := 43;
SiteLoc[47].y := 600;
SiteLoc[48].x := 43;
SiteLoc[48].y := 650;
SiteLoc[49].x := 43;
SiteLoc[49].y := 700;
SiteLoc[50].x := 43;
SiteLoc[50].y := 750;
SiteLoc[51].x := 43;
SiteLoc[51].y := 850;
SiteLoc[52].x := 43;
SiteLoc[52].y := 900;
SiteLoc[53].x := 43;
SiteLoc[53].y := 800;
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


#SURF ElectrodeTypes.pas definition for 2b polytrode
'''
Name := 'uMap54_2b';
Description := 'uMap54_2b, 50um spacing';
SiteSize.x := 15;
SiteSize.y := 15;
RoundSite := TRUE;
Created := FALSE;

NumPoints := 8;
Outline[0].x := -100;
Outline[0].y := -50;
Outline[1].x := -100;
Outline[1].y := 1350;
Outline[2].x := -20;
Outline[2].y := 1750;
Outline[3].x := 0;
Outline[3].y := 1850;
Outline[4].x := 20;
Outline[4].y := 1750;
Outline[5].x := 100;
Outline[5].y := 1350;
Outline[6].x := 100;
Outline[6].y := -50;
Outline[7].x := Outline[0].x;
Outline[7].y := Outline[0].y;

NumSites := 54;
CenterX := 0;
SiteLoc[0].x := -25;
SiteLoc[0].y := 1275;
SiteLoc[1].x := -25;
SiteLoc[1].y := 1175;
SiteLoc[2].x := -25;
SiteLoc[2].y := 1075;
SiteLoc[3].x := -25;
SiteLoc[3].y := 975;
SiteLoc[4].x := -25;
SiteLoc[4].y := 875;
SiteLoc[5].x := -25;
SiteLoc[5].y := 775;
SiteLoc[6].x := -25;
SiteLoc[6].y := 725;
SiteLoc[7].x := -25;
SiteLoc[7].y := 675;
SiteLoc[8].x := -25;
SiteLoc[8].y := 625;
SiteLoc[9].x := -25;
SiteLoc[9].y := 575;
SiteLoc[10].x := -25;
SiteLoc[10].y := 525;
SiteLoc[11].x := -25;
SiteLoc[11].y := 475;
SiteLoc[12].x := -25;
SiteLoc[12].y := 425;
SiteLoc[13].x := -25;
SiteLoc[13].y := 375;
SiteLoc[14].x := -25;
SiteLoc[14].y := 325;
SiteLoc[15].x := -25;
SiteLoc[15].y := 275;
SiteLoc[16].x := -25;
SiteLoc[16].y := 225;
SiteLoc[17].x := -25;
SiteLoc[17].y := 175;
SiteLoc[18].x := -25;
SiteLoc[18].y := 125;
SiteLoc[19].x := -25;
SiteLoc[19].y := 75;
SiteLoc[20].x := -25;
SiteLoc[20].y := 25;
SiteLoc[21].x := -25;
SiteLoc[21].y := 825;
SiteLoc[22].x := -25;
SiteLoc[22].y := 925;
SiteLoc[23].x := -25;
SiteLoc[23].y := 1025;
SiteLoc[24].x := -25;
SiteLoc[24].y := 1225;
SiteLoc[25].x := -25;
SiteLoc[25].y := 1325;
SiteLoc[26].x := -25;
SiteLoc[26].y := 1125;
SiteLoc[27].x := 25;
SiteLoc[27].y := 1300;
SiteLoc[28].x := 25;
SiteLoc[28].y := 1200;
SiteLoc[29].x := 25;
SiteLoc[29].y := 1100;
SiteLoc[30].x := 25;
SiteLoc[30].y := 1000;
SiteLoc[31].x := 25;
SiteLoc[31].y := 900;
SiteLoc[32].x := 25;
SiteLoc[32].y := 0;
SiteLoc[33].x := 25;
SiteLoc[33].y := 50;
SiteLoc[34].x := 25;
SiteLoc[34].y := 100;
SiteLoc[35].x := 25;
SiteLoc[35].y := 150;
SiteLoc[36].x := 25;
SiteLoc[36].y := 200;
SiteLoc[37].x := 25;
SiteLoc[37].y := 250;
SiteLoc[38].x := 25;
SiteLoc[38].y := 300;
SiteLoc[39].x := 25;
SiteLoc[39].y := 350;
SiteLoc[40].x := 25;
SiteLoc[40].y := 400;
SiteLoc[41].x := 25;
SiteLoc[41].y := 450;
SiteLoc[42].x := 25;
SiteLoc[42].y := 500;
SiteLoc[43].x := 25;
SiteLoc[43].y := 550;
SiteLoc[44].x := 25;
SiteLoc[44].y := 600;
SiteLoc[45].x := 25;
SiteLoc[45].y := 650;
SiteLoc[46].x := 25;
SiteLoc[46].y := 700;
SiteLoc[47].x := 25;
SiteLoc[47].y := 750;
SiteLoc[48].x := 25;
SiteLoc[48].y := 800;
SiteLoc[49].x := 25;
SiteLoc[49].y := 850;
SiteLoc[50].x := 25;
SiteLoc[50].y := 950;
SiteLoc[51].x := 25;
SiteLoc[51].y := 1150;
SiteLoc[52].x := 25;
SiteLoc[52].y := 1250;
SiteLoc[53].x := 25;
SiteLoc[53].y := 1050;

'''
