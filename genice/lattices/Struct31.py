# coding: utf-8
"""
Data source: Dutour Sikirić, Mathieu, Olaf Delgado-Friedrichs, and Michel Deza. “Space Fullerenes: a Computer Search for New Frank-Kasper Structures” Acta Crystallographica Section A Foundations of Crystallography 66.Pt 5 (2010): 602–615.

Cage composition:
 (12,14,15,16) = (20,8,8,4,)
"""

import genice.lattices
from genice.cell import cellvectors

class Lattice(genice.lattices.Lattice):
    def __init__(self):
        self.pairs="""
        1 225
        53 222
        32 83
        125 188
        66 54
        96 40
        158 94
        125 104
        223 12
        188 168
        13 206
        69 30
        96 57
        154 141
        158 108
        207 177
        165 149
        97 112
        97 110
        173 216
        134 19
        155 152
        148 78
        115 222
        66 212
        187 157
        54 2
        100 64
        197 131
        68 107
        80 152
        6 58
        69 111
        170 151
        140 136
        56 48
        148 82
        92 114
        21 159
        98 80
        148 95
        184 56
        145 85
        177 9
        111 58
        27 84
        180 133
        22 110
        46 223
        168 176
        100 23
        10 155
        112 210
        88 189
        22 32
        226 175
        182 221
        153 145
        183 223
        190 99
        61 12
        18 90
        212 117
        224 96
        213 117
        65 48
        61 100
        141 113
        144 115
        14 161
        62 221
        130 192
        129 67
        108 0
        88 90
        204 133
        4 109
        76 37
        122 105
        5 24
        174 97
        22 179
        226 157
        15 56
        86 12
        72 16
        120 144
        199 63
        199 62
        198 188
        88 84
        6 26
        21 213
        87 206
        182 183
        25 24
        208 211
        81 82
        68 206
        111 138
        83 102
        28 71
        218 103
        141 144
        32 138
        20 66
        34 187
        200 23
        142 223
        1 107
        158 207
        225 2
        107 120
        74 128
        162 0
        209 39
        173 19
        26 160
        33 104
        27 78
        203 45
        191 80
        31 43
        155 161
        153 131
        209 49
        99 75
        114 117
        17 95
        189 45
        193 156
        67 224
        60 126
        112 179
        91 59
        180 150
        116 30
        111 95
        153 186
        29 65
        28 65
        43 186
        210 166
        52 45
        52 42
        176 81
        34 82
        197 214
        59 204
        10 172
        148 210
        134 191
        103 55
        143 194
        158 211
        196 100
        196 99
        57 25
        86 217
        46 104
        27 45
        226 81
        219 169
        122 15
        35 60
        164 47
        26 77
        51 36
        162 207
        119 133
        11 61
        187 198
        46 90
        50 204
        118 58
        195 139
        191 36
        83 149
        164 9
        219 211
        108 132
        221 150
        79 152
        14 167
        4 172
        142 203
        181 30
        174 157
        142 11
        101 213
        153 41
        10 98
        108 169
        44 23
        35 212
        51 93
        178 0
        121 75
        208 215
        43 160
        8 221
        169 71
        140 224
        147 51
        101 70
        77 149
        117 216
        147 218
        75 62
        128 109
        17 174
        58 210
        127 159
        49 8
        201 55
        171 50
        181 102
        152 115
        84 198
        5 0
        27 168
        5 123
        151 35
        75 183
        8 128
        130 15
        202 93
        71 178
        124 7
        7 42
        129 160
        78 135
        201 163
        44 12
        33 121
        165 214
        18 196
        118 102
        109 39
        170 70
        89 97
        107 3
        77 41
        147 126
        91 62
        227 139
        162 139
        189 220
        220 188
        69 179
        142 121
        146 141
        178 143
        3 16
        208 137
        37 48
        33 199
        176 135
        101 127
        173 185
        59 8
        205 136
        146 161
        123 193
        21 92
        84 82
        31 57
        67 85
        21 130
        211 143
        154 201
        7 69
        220 89
        86 99
        168 175
        147 163
        73 139
        143 56
        123 132
        73 215
        167 222
        50 53
        29 184
        227 94
        22 6
        34 17
        208 37
        18 33
        17 7
        74 72
        119 63
        184 105
        78 42
        192 169
        219 76
        146 3
        5 194
        190 133
        31 131
        174 166
        119 182
        151 213
        125 217
        64 63
        118 77
        52 34
        122 212
        98 3
        123 215
        26 30
        191 202
        205 193
        1 79
        118 186
        200 109
        102 179
        13 113
        85 181
        187 203
        81 166
        93 54
        14 144
        15 71
        119 128
        64 183
        28 35
        140 164
        86 90
        92 76
        173 36
        73 48
        74 50
        89 176
        105 114
        31 136
        9 24
        202 225
        177 193
        120 80
        40 162
        32 160
        124 135
        216 55
        216 54
        47 41
        205 214
        145 214
        113 98
        106 120
        13 115
        70 163
        51 20
        218 2
        116 138
        200 63
        171 38
        209 172
        186 138
        68 72
        20 185
        134 79
        106 103
        137 132
        124 166
        23 150
        1 13
        215 207
        52 89
        218 127
        122 159
        209 53
        205 25
        92 60
        19 163
        49 72
        88 226
        10 206
        29 192
        85 136
        167 38
        224 165
        150 39
        113 202
        44 190
        170 105
        204 39
        192 137
        121 217
        140 156
        185 114
        197 177
        132 227
        74 87
        101 55
        87 4
        155 53
        44 182
        165 47
        38 59
        185 70
        19 2
        28 76
        154 106
        103 36
        170 126
        161 171
        172 16
        67 83
        110 95
        137 195
        189 217
        129 57
        18 11
        129 41
        73 194
        87 222
        219 184
        125 61
        116 149
        131 164
        4 38
        43 181
        94 9
        66 127
        167 49
        116 145
        178 195
        14 68
        154 134
        91 190
        40 24
        180 199
        46 175
        159 126
        180 196
        42 112
        29 151
        96 197
        156 40
        106 225
        201 93
        6 124
        20 60
        171 16
        135 110
        200 91
        25 47
        130 37
        104 64
        94 194
        203 175
        11 198
        146 79
        65 195
        220 157
        156 227
        """

        self.waters="""
        0.45717 0.91434 0.60814
        0.33323 0.04161 0.30662
        0.33566 0.04283 0.37
        0.0 0.0 0.27211
        0.70839 0.04162 0.19338
        0.29283 0.95717 0.63
        0.41667 0.08333 0.83856
        0.04283 0.70717 0.87
        0.08333 0.66667 0.16145
        0.41667 0.33333 0.66145
        0.58345 0.79173 0.25
        0.70717 0.66434 0.01853
        0.0 0.0 0.05559
        0.54162 0.08323 0.28451
        0.08344 0.54172 0.25
        0.29283 0.95717 0.51853
        0.0 0.0 0.22789
        0.87384 0.74767 0.89199
        0.54283 0.45717 0.03706
        0.33566 0.29283 0.37
        0.581 0.7905 0.42802
        0.95717 0.66434 0.48147
        0.66667 0.08333 0.83856
        0.74767 0.87384 0.10802
        0.41667 0.08333 0.66145
        0.29162 0.95839 0.69338
        0.29161 0.95839 0.80662
        0.33334 0.66667 0.94877
        0.66434 0.95717 0.51853
        0.04283 0.33566 0.51853
        0.04162 0.70839 0.80662
        0.33333 0.66667 0.72802
        0.66677 0.95839 0.80662
        0.66667 0.33334 0.05123
        0.79051 0.58101 0.92802
        0.70717 0.04283 0.48147
        0.70717 0.66434 0.37
        0.87384 0.74767 0.54142
        0.95839 0.29162 0.19338
        0.58333 0.66667 0.16145
        0.66667 0.08333 0.66145
        0.41657 0.20828 0.75
        0.08566 0.54283 0.89186
        0.33333 0.66667 0.77198
        0.0 0.0 0.09598
        0.08567 0.54284 0.96294
        0.33566 0.04283 0.01853
        0.2499 0.12495 0.71536
        0.0 0.0 0.55559
        0.95838 0.66677 0.19338
        0.33323 0.29162 0.19338
        0.74767 0.87384 0.39199
        0.919 0.4595 0.92802
        0.54161 0.45839 0.21549
        0.54283 0.08566 0.39186
        0.91434 0.45717 0.39186
        0.25233 0.12616 0.54142
        0.45839 0.91678 0.71549
        0.41667 0.33333 0.83856
        0.08333 0.41667 0.16145
        0.74767 0.87384 0.45858
        0.74768 0.87384 0.04142
        0.91434 0.45717 0.10814
        0.54283 0.08566 0.10814
        0.4595 0.919 0.07199
        0.87384 0.12617 0.54142
        0.4595 0.919 0.42802
        0.79172 0.20828 0.75
        0.20828 0.79172 0.25
        0.91667 0.58333 0.83856
        0.2095 0.41899 0.42802
        0.45717 0.91434 0.53706
        0.12494 0.87505 0.21536
        0.0 0.0 0.59598
        0.33323 0.04161 0.19338
        0.08101 0.54051 0.07199
        0.66434 0.70717 0.51853
        0.24989 0.12495 0.78464
        0.33333 0.66667 0.90389
        0.33323 0.29161 0.30662
        0.70839 0.66678 0.30662
        0.41899 0.2095 0.92802
        0.54051 0.4595 0.92802
        0.87505 0.12495 0.78464
        0.45717 0.54283 0.96294
        0.91656 0.45828 0.75
        0.12617 0.25233 0.04142
        0.54162 0.08323 0.21549
        0.29283 0.33566 0.98147
        0.7905 0.2095 0.92802
        0.33567 0.29283 0.01853
        0.95717 0.29283 0.13
        0.70717 0.66434 0.48147
        0.70717 0.04283 0.37
        0.29283 0.33566 0.63
        0.66434 0.70717 0.87
        0.66677 0.95839 0.69338
        0.87383 0.12617 0.89199
        0.75011 0.87505 0.28464
        0.2095 0.419 0.07199
        0.581 0.7905 0.07199
        0.08101 0.54051 0.42802
        0.04161 0.33323 0.80662
        0.95717 0.66434 0.37
        0.54284 0.08567 0.03706
        0.33566 0.29283 0.48147
        0.08333 0.66667 0.33856
        0.12495 0.87505 0.28464
        0.33333 0.66666 0.59611
        0.58333 0.91667 0.16145
        0.66434 0.95717 0.87
        0.66667 0.58333 0.83856
        0.04283 0.33566 0.87
        0.70839 0.04161 0.30662
        0.54283 0.45717 0.46294
        0.66667 0.33334 0.27198
        0.87505 0.75011 0.78464
        0.66667 0.33333 0.44877
        0.29161 0.33323 0.80662
        0.33566 0.04283 0.13
        0.95839 0.66678 0.30662
        0.91434 0.45717 0.03706
        0.33566 0.04283 0.48147
        0.04283 0.70717 0.63
        0.29283 0.95717 0.87
        0.70717 0.04284 0.01853
        0.0 0.0 0.44441
        0.2095 0.79051 0.42802
        0.33333 0.91667 0.16145
        0.54172 0.08344 0.75
        0.04283 0.70717 0.51853
        0.45839 0.54161 0.71549
        0.08566 0.54283 0.60814
        0.33566 0.29283 0.13
        0.33333 0.41667 0.33856
        0.45717 0.91434 0.89186
        0.08323 0.54161 0.71549
        0.919 0.4595 0.57199
        0.66678 0.70839 0.80662
        0.87384 0.12617 0.60802
        0.04161 0.33323 0.69338
        0.95839 0.29161 0.30662
        0.95717 0.66434 0.01853
        0.41899 0.20949 0.57199
        0.91677 0.45839 0.28451
        0.79172 0.58344 0.75
        0.12495 0.24989 0.28464
        0.0 0.0 0.40402
        0.45717 0.54283 0.89186
        0.0 0.0 0.77211
        0.70717 0.66434 0.13
        0.95717 0.29283 0.48147
        0.54161 0.45839 0.28451
        0.54172 0.45828 0.75
        0.08333 0.41667 0.33856
        0.45828 0.54172 0.25
        0.91667 0.33333 0.66145
        0.0 0.0 0.94441
        0.45717 0.54283 0.60814
        0.12617 0.87384 0.45858
        0.45839 0.91677 0.78451
        0.20828 0.41656 0.25
        0.66434 0.95717 0.63
        0.12617 0.25233 0.39199
        0.29162 0.33323 0.69338
        0.0 0.0 0.72789
        0.25233 0.12617 0.89199
        0.91678 0.45839 0.21549
        0.45717 0.91434 0.96294
        0.33333 0.66667 0.55123
        0.12617 0.25233 0.45858
        0.12495 0.24989 0.21536
        0.75011 0.87506 0.21536
        0.54283 0.45717 0.39186
        0.0 0.0 0.90402
        0.29283 0.95717 0.98147
        0.5405 0.08101 0.92802
        0.66667 0.58333 0.66145
        0.5405 0.081 0.57199
        0.91666 0.33333 0.83856
        0.54283 0.45717 0.10814
        0.08323 0.54162 0.78451
        0.12617 0.87384 0.10802
        0.2095 0.7905 0.07199
        0.29283 0.33566 0.51853
        0.4595 0.5405 0.42802
        0.45838 0.54161 0.78451
        0.87384 0.74768 0.95858
        0.66434 0.95717 0.98147
        0.04283 0.33566 0.98147
        0.12617 0.25233 0.10802
        0.58333 0.66667 0.33856
        0.08566 0.54283 0.53706
        0.91667 0.58333 0.66145
        0.25233 0.12616 0.60802
        0.7905 0.2095 0.57199
        0.4595 0.54051 0.07199
        0.66678 0.70839 0.69338
        0.66434 0.70717 0.98147
        0.66667 0.33333 0.09611
        0.70717 0.04283 0.13
        0.95717 0.29283 0.37
        0.58333 0.91667 0.33856
        0.04284 0.70717 0.98147
        0.33333 0.41667 0.16145
        0.04161 0.70838 0.69338
        0.45829 0.91656 0.25
        0.66434 0.70717 0.63
        0.7905 0.58101 0.57199
        0.70839 0.66677 0.19338
        0.29283 0.33566 0.87
        0.54051 0.4595 0.57199
        0.54283 0.08566 0.46294
        0.91434 0.45717 0.46294
        0.87505 0.7501 0.71536
        0.87384 0.74767 0.60802
        0.66667 0.33333 0.40389
        0.95717 0.29283 0.01853
        0.12617 0.87384 0.39199
        0.45717 0.54283 0.53706
        0.87383 0.12617 0.95858
        0.95717 0.66434 0.13
        0.66667 0.33334 0.22802
        0.12616 0.87383 0.04142
        0.87505 0.12495 0.71536
        0.33333 0.91667 0.33856
        0.25233 0.12617 0.95858
        0.04283 0.33566 0.63
        """

        self.coord= "relative"

        self.cages="""
        12 0.34265 0.17132 -0.57412
        14 0.33333 0.66667 0.12226
        16 -0.66667 -0.33333 -0.51744
        14 -0.33333 -0.66667 0.69989
        14 -0.66667 -0.33333 -0.62226
        15 0.0 1.0 -0.66156
        12 0.17132 -0.17133 0.57412
        12 -0.16645 0.16646 0.25
        12 0.16645 0.33291 -0.25
        15 -0.66667 -0.33333 0.66207
        12 -0.16646 -0.33291 0.25
        12 -0.17133 0.17132 -0.57412
        14 0.33333 0.66667 0.19989
        12 0.16645 -0.16646 0.75
        15 0.0 -1.0 -1.16156
        15 0.33333 0.66667 -0.16207
        16 -0.33333 -0.66667 0.51744
        15 0.66667 0.33333 0.16207
        15 -1.0 0.0 1.16156
        14 0.66667 0.33333 -0.19989
        12 -1.0 0.0 0.5
        12 -0.17132 0.17133 0.07412
        12 0.33291 0.16645 0.25
        15 -0.33333 -0.66667 -0.66207
        12 0.17133 0.34265 0.57412
        12 0.17133 -0.17132 -0.07412
        15 1.0 0.0 0.66156
        14 -0.33333 -0.66667 0.62226
        14 -0.66667 -0.33333 -0.69989
        12 1.0 0.0 0.0
        12 -0.17132 -0.34265 -0.57412
        14 0.66667 0.33333 -0.12226
        12 -0.34265 -0.17132 0.57412
        12 -0.33291 -0.16645 0.75
        12 0.34265 0.17132 0.07412
        12 -0.17133 -0.34265 0.07412
        12 0.17132 0.34265 -0.07412
        16 0.66667 0.33333 -0.01744
        12 -0.34265 -0.17132 -0.07412
        16 0.33333 0.66667 0.01744
        """

        self.bondlen = 3


        self.cell = """
        12.748400612611047 0.0 0.0
        -6.374200306305521 11.040438788142268 0.0
        4.488261564706488e-15 7.773897067730223e-15 73.29887389296881
        """

        self.density = 0.6605828002309204



        self.cell = cellvectors(a=12.748400612611047,
                           b=12.748400612611047,
                           c=73.29887389296881,
                           C=119.99999999999999)
