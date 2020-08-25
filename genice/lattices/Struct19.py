# coding: utf-8
"""
Data source: Dutour Sikirić, Mathieu, Olaf Delgado-Friedrichs, and Michel Deza. “Space Fullerenes: a Computer Search for New Frank-Kasper Structures” Acta Crystallographica Section A Foundations of Crystallography 66.Pt 5 (2010): 602–615.

Cage composition:
 (12,14,15,16) = (16,0,0,8,)
"""

import genice.lattices
from genice.cell import cellvectors

class Lattice(genice.lattices.Lattice):
    def __init__(self):
        self.pairs="""
        64 89
        67 79
        8 54
        18 75
        19 32
        124 111
        67 113
        122 58
        52 96
        34 63
        28 123
        110 133
        90 89
        100 78
        33 122
        135 21
        75 131
        41 66
        67 65
        46 64
        70 133
        82 114
        16 49
        84 115
        124 96
        27 9
        77 97
        134 117
        66 116
        77 116
        81 50
        27 89
        13 85
        25 95
        31 94
        18 134
        121 72
        92 113
        102 126
        111 129
        1 34
        16 123
        32 109
        104 76
        103 4
        104 6
        105 6
        81 98
        37 78
        125 128
        73 128
        70 127
        1 59
        74 128
        83 102
        91 103
        54 106
        93 59
        62 15
        3 5
        108 112
        47 109
        57 111
        17 58
        85 126
        113 129
        13 2
        35 121
        3 17
        92 41
        56 61
        55 60
        120 127
        120 126
        37 131
        55 27
        53 58
        94 118
        52 26
        53 59
        122 97
        104 64
        19 88
        28 55
        47 21
        112 22
        30 80
        63 71
        62 73
        7 130
        3 132
        40 134
        83 48
        57 50
        91 51
        27 75
        25 68
        62 84
        3 97
        64 83
        81 99
        40 98
        41 100
        9 115
        94 86
        96 87
        97 88
        12 13
        28 29
        32 77
        134 112
        124 10
        47 82
        20 109
        11 92
        70 118
        100 101
        37 74
        107 119
        93 108
        24 110
        34 4
        109 78
        110 80
        108 80
        55 46
        76 61
        32 61
        24 127
        12 105
        35 84
        36 85
        11 103
        2 72
        75 10
        41 79
        38 66
        133 23
        39 113
        14 135
        11 70
        0 83
        119 102
        2 84
        14 69
        12 73
        0 82
        93 90
        18 106
        7 112
        52 129
        37 62
        110 132
        17 135
        33 121
        25 4
        14 15
        74 101
        56 60
        13 49
        107 117
        92 51
        0 25
        16 65
        9 30
        30 2
        132 107
        71 99
        130 57
        77 100
        44 79
        44 78
        45 80
        45 72
        99 23
        63 50
        17 82
        31 98
        133 66
        33 101
        123 79
        125 129
        38 60
        31 117
        35 58
        101 51
        103 69
        36 1
        106 74
        20 124
        44 21
        45 23
        45 22
        44 15
        36 90
        38 54
        69 127
        7 59
        46 61
        52 73
        14 35
        7 122
        90 102
        0 5
        71 72
        130 88
        19 39
        108 107
        40 10
        18 60
        71 85
        38 123
        19 111
        131 115
        118 126
        132 116
        94 119
        34 130
        104 65
        11 86
        12 87
        54 23
        98 6
        120 48
        53 114
        120 49
        56 117
        56 116
        53 115
        125 50
        125 51
        65 48
        20 114
        31 76
        63 91
        1 95
        26 68
        39 86
        42 5
        42 4
        36 87
        43 6
        20 131
        43 10
        106 22
        68 48
        99 118
        46 47
        26 15
        28 21
        33 22
        68 69
        95 96
        29 30
        9 93
        119 5
        40 57
        26 67
        24 29
        8 81
        95 114
        29 49
        121 91
        8 105
        39 76
        24 135
        42 86
        42 88
        43 89
        43 87
        8 128
        16 105
        """

        self.waters="""
        0.5625 0.1875 0.375
        0.625 0.8125 0.4375
        0.0625 0.6875 0.375
        0.375 0.1875 0.5625
        0.71875 0.03125 0.53125
        0.53125 0.21875 0.53125
        0.71875 0.53125 0.03125
        0.4375 0.8125 0.625
        0.9375 0.625 0.9375
        0.3125 0.625 0.3125
        0.53125 0.71875 0.03125
        0.875 0.1875 0.6875
        0.875 0.6875 0.1875
        0.9375 0.625 0.3125
        0.0625 0.0625 0.375
        0.03125 0.03125 0.21875
        0.9375 0.4375 0.125
        0.3125 0.125 0.4375
        0.3125 0.625 0.9375
        0.5625 0.0625 0.875
        0.4375 0.9375 0.125
        0.21875 0.21875 0.21875
        0.21875 0.71875 0.71875
        0.03125 0.53125 0.71875
        0.125 0.3125 0.4375
        0.6875 0.0625 0.375
        0.875 0.0625 0.1875
        0.375 0.5625 0.1875
        0.1875 0.375 0.1875
        0.125 0.4375 0.3125
        0.1875 0.5625 0.375
        0.625 0.4375 0.8125
        0.4375 0.125 0.9375
        0.1875 0.875 0.6875
        0.6875 0.875 0.5625
        0.125 0.9375 0.4375
        0.6875 0.6875 0.375
        0.1875 0.875 0.0625
        0.125 0.4375 0.9375
        0.6875 0.1875 0.875
        0.5625 0.6875 0.875
        0.0625 0.1875 0.875
        0.625 0.125 0.625
        0.625 0.625 0.125
        0.125 0.125 0.125
        0.125 0.625 0.625
        0.4375 0.3125 0.125
        0.375 0.1875 0.1875
        0.78125 0.28125 0.28125
        0.96875 0.46875 0.28125
        0.78125 0.78125 0.78125
        0.96875 0.96875 0.78125
        0.8125 0.9375 0.125
        0.375 0.875 0.375
        0.0625 0.5625 0.875
        0.3125 0.4375 0.125
        0.375 0.375 0.875
        0.625 0.8125 0.8125
        0.28125 0.96875 0.46875
        0.46875 0.78125 0.46875
        0.28125 0.46875 0.96875
        0.46875 0.28125 0.96875
        0.0625 0.875 0.1875
        0.8125 0.8125 0.625
        0.5625 0.375 0.1875
        0.8125 0.3125 0.125
        0.125 0.3125 0.8125
        0.875 0.1875 0.0625
        0.8125 0.125 0.3125
        0.9375 0.125 0.4375
        0.9375 0.3125 0.625
        0.875 0.6875 0.5625
        0.03125 0.71875 0.53125
        0.9375 0.8125 0.125
        0.125 0.8125 0.9375
        0.375 0.6875 0.0625
        0.625 0.3125 0.9375
        0.3125 0.125 0.8125
        0.21875 0.03125 0.03125
        0.03125 0.21875 0.03125
        0.21875 0.53125 0.53125
        0.8125 0.625 0.8125
        0.4375 0.125 0.3125
        0.625 0.3125 0.3125
        0.125 0.8125 0.3125
        0.8125 0.625 0.4375
        0.71875 0.21875 0.71875
        0.71875 0.71875 0.21875
        0.53125 0.03125 0.71875
        0.53125 0.53125 0.21875
        0.5625 0.5625 0.375
        0.9375 0.9375 0.625
        0.9375 0.125 0.8125
        0.4375 0.625 0.4375
        0.6875 0.375 0.6875
        0.625 0.9375 0.3125
        0.6875 0.875 0.1875
        0.375 0.0625 0.6875
        0.6875 0.5625 0.875
        0.875 0.5625 0.6875
        0.1875 0.0625 0.875
        0.125 0.9375 0.8125
        0.625 0.4375 0.4375
        0.875 0.0625 0.5625
        0.6875 0.375 0.0625
        0.875 0.5625 0.0625
        0.1875 0.6875 0.875
        0.4375 0.4375 0.625
        0.375 0.5625 0.5625
        0.375 0.0625 0.0625
        0.1875 0.375 0.5625
        0.625 0.9375 0.9375
        0.375 0.6875 0.6875
        0.8125 0.125 0.9375
        0.46875 0.96875 0.28125
        0.28125 0.78125 0.28125
        0.28125 0.28125 0.78125
        0.46875 0.46875 0.78125
        0.8125 0.4375 0.625
        0.5625 0.375 0.5625
        0.875 0.375 0.375
        0.0625 0.875 0.5625
        0.3125 0.9375 0.625
        0.0625 0.375 0.0625
        0.5625 0.875 0.0625
        0.875 0.875 0.875
        0.78125 0.46875 0.46875
        0.96875 0.28125 0.46875
        0.96875 0.78125 0.96875
        0.78125 0.96875 0.96875
        0.5625 0.875 0.6875
        0.3125 0.8125 0.125
        0.3125 0.3125 0.625
        0.0625 0.375 0.6875
        0.4375 0.625 0.8125
        0.1875 0.1875 0.375
        """

        self.coord= "relative"

        self.cages="""
        12 0.75 0.25 0.5
        16 -0.125 -0.125 -0.625
        12 0.5 0.0 0.5
        12 -0.25 -0.5 -0.75
        16 0.125 0.625 0.125
        12 0.5 -0.75 -0.25
        12 -0.5 0.75 0.25
        16 -0.375 0.625 -0.375
        12 1.25 0.25 1.0
        12 0.25 0.5 0.75
        16 0.625 0.125 1.125
        12 0.0 0.5 0.5
        16 1.375 0.375 0.375
        12 0.75 0.75 1.0
        12 1.0 0.25 1.25
        16 0.875 0.375 -0.125
        12 -0.75 -0.25 -0.5
        12 0.0 0.75 0.75
        16 -0.625 -0.125 -1.125
        12 0.25 1.0 1.25
        16 0.125 0.125 0.625
        12 -0.25 0.0 -0.25
        12 0.5 0.5 0.0
        12 0.0 0.0 0.0
        """

        self.bondlen = 3


        self.cell = """
        18.028244506884192 18.028244506884192 18.028244506884192
        """

        self.density = 0.6937617372557815



        self.cell = cellvectors(a=18.028244506884192,
                           b=18.028244506884192,
                           c=18.028244506884192)
