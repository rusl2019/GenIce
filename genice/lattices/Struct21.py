# coding: utf-8
"""
Data source: Dutour Sikirić, Mathieu, Olaf Delgado-Friedrichs, and Michel Deza. “Space Fullerenes: a Computer Search for New Frank-Kasper Structures” Acta Crystallographica Section A Foundations of Crystallography 66.Pt 5 (2010): 602–615.

Cage composition:
 (12,14,15,16) = (16,20,8,0,)
"""

import genice.lattices
from genice.cell import cellvectors

class Lattice(genice.lattices.Lattice):
    def __init__(self):
        self.pairs="""
        155 80
        106 222
        250 197
        220 159
        137 101
        106 97
        34 25
        142 184
        93 50
        212 234
        251 150
        49 123
        31 18
        184 177
        177 183
        27 32
        43 79
        155 143
        74 64
        239 250
        100 39
        208 217
        71 240
        182 152
        193 116
        133 91
        89 41
        24 114
        34 203
        163 162
        26 164
        120 52
        3 161
        37 206
        19 51
        204 59
        101 33
        75 84
        16 94
        62 170
        118 114
        221 218
        183 6
        32 222
        104 160
        42 245
        6 238
        21 104
        43 48
        161 92
        169 175
        5 88
        87 219
        37 201
        49 73
        71 85
        29 138
        54 38
        129 217
        17 65
        227 128
        149 156
        63 117
        96 144
        1 139
        221 94
        88 20
        32 233
        198 194
        37 90
        55 83
        89 244
        241 74
        25 104
        73 124
        139 105
        236 105
        83 7
        244 165
        9 61
        94 21
        199 135
        237 150
        149 242
        112 13
        89 230
        55 220
        4 44
        212 19
        243 181
        77 214
        237 113
        2 97
        81 210
        15 143
        28 12
        3 61
        23 176
        178 215
        22 245
        216 143
        16 191
        168 117
        22 70
        42 212
        147 251
        182 134
        147 191
        239 57
        139 190
        168 153
        5 66
        197 52
        202 39
        66 116
        184 35
        125 170
        201 207
        18 121
        4 179
        121 108
        244 223
        69 162
        0 63
        136 145
        216 41
        198 82
        80 223
        98 44
        20 223
        189 87
        56 139
        133 34
        117 156
        205 161
        232 25
        35 64
        43 133
        68 219
        169 75
        70 19
        207 148
        70 214
        218 72
        1 37
        209 84
        47 40
        93 234
        216 165
        48 21
        77 138
        10 207
        118 187
        132 151
        2 248
        107 237
        136 30
        209 248
        27 39
        167 51
        49 176
        1 227
        209 98
        47 36
        87 196
        10 62
        162 134
        42 110
        81 241
        76 183
        65 190
        21 205
        224 164
        31 166
        185 250
        22 120
        126 79
        229 17
        11 98
        132 165
        236 247
        228 145
        28 69
        25 88
        127 138
        122 144
        123 145
        214 57
        55 145
        54 107
        152 35
        153 36
        27 171
        178 78
        155 48
        29 96
        225 59
        133 225
        232 244
        125 121
        16 9
        213 109
        141 180
        208 45
        136 8
        75 195
        218 111
        193 227
        231 122
        16 160
        245 186
        78 77
        153 33
        15 45
        198 211
        189 85
        0 38
        208 66
        45 180
        225 126
        28 231
        95 199
        10 190
        114 175
        164 241
        142 210
        52 246
        188 140
        103 233
        107 33
        220 158
        210 224
        134 101
        97 105
        142 182
        212 235
        251 149
        15 230
        191 40
        192 38
        130 222
        119 58
        23 12
        31 67
        163 137
        102 6
        201 65
        62 108
        31 27
        23 102
        235 85
        84 234
        76 82
        8 57
        58 7
        180 148
        91 165
        206 99
        40 211
        176 146
        121 65
        19 50
        203 59
        141 91
        249 236
        60 237
        43 217
        38 113
        151 204
        125 95
        167 235
        48 111
        160 92
        236 233
        172 200
        220 239
        126 92
        127 158
        157 146
        32 90
        183 124
        54 40
        60 63
        166 249
        235 179
        146 96
        9 211
        171 17
        202 68
        128 230
        134 6
        90 17
        178 120
        249 190
        238 74
        232 151
        60 82
        47 194
        226 138
        174 204
        13 158
        12 159
        86 153
        26 238
        86 137
        221 156
        72 9
        49 164
        8 83
        166 108
        127 58
        142 69
        178 115
        130 125
        181 122
        228 159
        192 14
        0 194
        243 102
        154 230
        120 7
        81 36
        175 196
        50 186
        187 51
        186 52
        51 188
        99 116
        80 205
        5 89
        179 110
        157 73
        86 210
        226 159
        28 13
        172 115
        97 229
        131 4
        131 140
        226 30
        127 200
        185 187
        50 53
        112 146
        79 104
        155 132
        63 14
        177 194
        29 115
        242 211
        78 140
        72 149
        176 228
        106 56
        61 113
        103 106
        147 192
        157 163
        231 224
        58 239
        88 91
        189 2
        53 197
        130 193
        166 103
        205 232
        141 217
        172 144
        196 110
        109 195
        170 135
        215 8
        184 238
        167 114
        93 24
        243 46
        119 77
        174 61
        26 69
        219 229
        112 162
        169 240
        172 7
        103 248
        101 74
        11 179
        45 99
        181 13
        248 67
        68 247
        161 204
        203 160
        36 64
        219 67
        243 73
        83 197
        215 70
        167 131
        191 156
        213 24
        3 192
        234 44
        119 185
        147 92
        18 222
        124 137
        2 247
        111 151
        60 152
        129 199
        216 180
        107 242
        173 213
        175 85
        174 218
        100 195
        72 203
        136 115
        201 95
        118 110
        87 195
        24 75
        1 154
        42 109
        4 173
        76 182
        241 124
        34 111
        169 39
        5 135
        71 247
        224 163
        102 112
        22 119
        188 250
        128 66
        126 80
        209 240
        96 158
        229 18
        0 35
        79 20
        100 67
        170 148
        46 12
        199 41
        206 148
        108 56
        100 98
        202 11
        173 246
        208 20
        44 109
        99 154
        47 117
        168 82
        193 56
        130 206
        207 41
        15 223
        55 144
        171 68
        157 122
        95 227
        181 200
        46 30
        233 240
        214 188
        213 186
        251 59
        221 14
        177 81
        90 105
        198 113
        129 128
        226 57
        249 171
        118 173
        132 225
        46 123
        202 196
        141 135
        62 116
        23 26
        78 246
        174 150
        14 150
        215 53
        231 123
        10 154
        86 76
        152 33
        245 187
        140 53
        189 84
        29 228
        54 64
        246 185
        131 93
        3 94
        129 143
        168 242
        71 11
        200 30
        """

        self.waters="""
        0.3125 0.05868 0.29674
        0.30868 0.0625 0.95326
        0.25 0.25 0.82733
        0.1875 0.05868 0.20326
        0.0 0.75 0.71034
        0.0 0.375 0.03966
        0.0 0.05868 0.39719
        0.68368 0.375 0.57799
        0.31632 0.375 0.57799
        0.0 0.375 0.22309
        0.0 0.875 0.96034
        0.0 0.875 0.77691
        0.31632 0.0625 0.48017
        0.68368 0.0625 0.48017
        0.5 0.05868 0.23659
        0.30868 0.9375 0.04674
        0.1875 0.44132 0.20326
        0.69132 0.06632 0.88577
        0.69132 0.43368 0.88577
        0.25 0.25 0.67267
        0.18368 0.625 0.07799
        0.30868 0.25 0.14719
        0.875 0.25 0.63708
        0.19132 0.9375 0.45326
        0.5 0.625 0.72309
        0.0 0.375 0.11293
        0.31632 0.875 0.42201
        0.69132 0.75 0.85282
        0.5 0.125 0.46034
        0.0 0.75 0.53966
        0.125 0.25 0.52691
        0.81632 0.55868 0.86424
        0.5 0.75 0.875
        0.75 0.75 0.32733
        0.81632 0.44132 0.13577
        0.375 0.93368 0.32799
        0.375 0.56632 0.32799
        0.5 0.9375 0.93986
        0.1875 0.93368 0.26983
        0.6875 0.75 0.81015
        0.1875 0.56632 0.26983
        0.81632 0.0625 0.01983
        0.875 0.25 0.71034
        0.5 0.44132 0.10282
        0.0 0.55868 0.73659
        0.375 0.75 0.02691
        0.19132 0.25 0.48659
        0.3125 0.44132 0.29674
        0.5 0.25 0.125
        0.19132 0.5625 0.45326
        0.375 0.43368 0.67201
        0.375 0.06632 0.67201
        0.68368 0.55868 0.63577
        0.31632 0.55868 0.63577
        0.125 0.75 0.28966
        0.5 0.625 0.53966
        0.18368 0.375 0.92201
        0.31632 0.125 0.57799
        0.68368 0.125 0.57799
        0.75 0.75 0.17267
        0.6875 0.05868 0.29674
        0.0 0.125 0.22309
        0.0 0.625 0.96034
        0.5 0.125 0.27691
        0.25 0.75 0.32733
        0.81632 0.125 0.92201
        0.18368 0.4375 0.01983
        0.875 0.43368 0.82799
        0.875 0.06632 0.82799
        0.5 0.0 0.42267
        0.125 0.25 0.63708
        0.1875 0.94132 0.79674
        0.8125 0.44132 0.20326
        0.0 0.4375 0.43986
        0.125 0.75 0.36293
        0.5 0.55868 0.76341
        0.80868 0.25 0.35282
        0.0 0.94132 0.60282
        0.0 0.75 0.625
        0.30868 0.56632 0.11424
        0.30868 0.93368 0.11424
        0.31632 0.44132 0.36424
        0.8125 0.25 0.31015
        0.5 0.5 0.57733
        0.3125 0.43368 0.76983
        0.3125 0.06632 0.76983
        0.68368 0.44132 0.36424
        0.625 0.25 0.78966
        0.0 0.5 0.07733
        0.0 0.125 0.03966
        0.5 0.94132 0.89719
        0.81632 0.625 0.07799
        0.25 0.75 0.17267
        0.3125 0.55868 0.70326
        0.3125 0.25 0.18986
        0.625 0.25 0.97309
        0.80868 0.75 0.51341
        0.375 0.25 0.86293
        0.0 0.625 0.77691
        0.30868 0.75 0.98659
        0.8125 0.55868 0.79674
        0.875 0.75 0.36293
        0.0 0.0625 0.43986
        0.18368 0.55868 0.86424
        0.18368 0.44132 0.13577
        0.30868 0.06632 0.88577
        0.30868 0.43368 0.88577
        0.875 0.75 0.28966
        0.0 0.5 0.92267
        0.8125 0.43368 0.73017
        0.8125 0.06632 0.73017
        0.69132 0.25 0.14719
        0.80868 0.9375 0.45326
        0.0 0.05868 0.26341
        0.5 0.875 0.72309
        0.0 0.5625 0.56015
        0.18368 0.5625 0.98017
        0.5 0.375 0.27691
        0.6875 0.94132 0.70326
        0.80868 0.06632 0.61424
        0.80868 0.43368 0.61424
        0.81632 0.375 0.92201
        0.68368 0.4375 0.48017
        0.31632 0.4375 0.48017
        0.0 0.44132 0.39719
        0.69132 0.4375 0.95326
        0.375 0.75 0.13708
        0.80868 0.0625 0.54674
        0.30868 0.25 0.01341
        0.5 0.25 0.03966
        0.5 0.5625 0.93986
        0.1875 0.75 0.68986
        0.69132 0.93368 0.11424
        0.69132 0.56632 0.11424
        0.80868 0.93368 0.38577
        0.81632 0.4375 0.01983
        0.19132 0.4375 0.54674
        0.80868 0.56632 0.38577
        0.0 0.9375 0.56015
        0.18368 0.125 0.92201
        0.19132 0.75 0.64719
        0.69132 0.5625 0.04674
        0.5 0.125 0.38708
        0.5 0.0625 0.06015
        0.68368 0.5625 0.51983
        0.31632 0.5625 0.51983
        0.875 0.75 0.47309
        0.375 0.75 0.21034
        0.69132 0.75 0.98659
        0.6875 0.56632 0.23017
        0.6875 0.93368 0.23017
        0.81632 0.05868 0.13577
        0.625 0.93368 0.32799
        0.625 0.56632 0.32799
        0.18368 0.9375 0.98017
        0.5 0.05868 0.10282
        0.5 0.44132 0.23659
        0.80868 0.5625 0.45326
        0.68368 0.9375 0.51983
        0.31632 0.9375 0.51983
        0.125 0.56632 0.17201
        0.125 0.93368 0.17201
        0.68368 0.875 0.42201
        0.68368 0.625 0.42201
        0.31632 0.625 0.42201
        0.81632 0.875 0.07799
        0.0 0.625 0.88708
        0.3125 0.94132 0.70326
        0.6875 0.44132 0.29674
        0.5 0.75 0.78966
        0.81632 0.5625 0.98017
        0.81632 0.94132 0.86424
        0.80868 0.4375 0.54674
        0.8125 0.75 0.68986
        0.8125 0.05868 0.20326
        0.5 0.94132 0.76341
        0.125 0.75 0.47309
        0.19132 0.25 0.35282
        0.0 0.55868 0.60282
        0.0 0.94132 0.73659
        0.625 0.75 0.02691
        0.80868 0.25 0.48659
        0.68368 0.05868 0.36424
        0.0 0.25 0.375
        0.31632 0.05868 0.36424
        0.68368 0.94132 0.63577
        0.625 0.43368 0.67201
        0.625 0.06632 0.67201
        0.31632 0.94132 0.63577
        0.375 0.25 0.78966
        0.0 0.0 0.92267
        0.3125 0.56632 0.23017
        0.3125 0.93368 0.23017
        0.30868 0.4375 0.95326
        0.1875 0.25 0.31015
        0.6875 0.43368 0.76983
        0.6875 0.06632 0.76983
        0.5 0.625 0.61293
        0.0 0.25 0.28966
        0.69132 0.25 0.01341
        0.875 0.25 0.52691
        0.69132 0.0625 0.95326
        0.8125 0.94132 0.79674
        0.875 0.56632 0.17201
        0.875 0.93368 0.17201
        0.18368 0.05868 0.13577
        0.5 0.75 0.96034
        0.81632 0.9375 0.98017
        0.30868 0.5625 0.04674
        0.1875 0.55868 0.79674
        0.5 0.375 0.38708
        0.0 0.44132 0.26341
        0.125 0.25 0.71034
        0.6875 0.55868 0.70326
        0.19132 0.06632 0.61424
        0.19132 0.43368 0.61424
        0.69132 0.9375 0.04674
        0.5 0.4375 0.06015
        0.6875 0.25 0.18986
        0.75 0.25 0.82733
        0.5 0.875 0.53966
        0.5 0.25 0.21034
        0.5 0.55868 0.89719
        0.18368 0.875 0.07799
        0.5 0.5 0.42267
        0.625 0.75 0.13708
        0.19132 0.0625 0.54674
        0.375 0.25 0.97309
        0.19132 0.75 0.51341
        0.625 0.25 0.86293
        0.18368 0.0625 0.01983
        0.5 0.375 0.46034
        0.0 0.125 0.11293
        0.30868 0.75 0.85282
        0.1875 0.43368 0.73017
        0.1875 0.06632 0.73017
        0.18368 0.94132 0.86424
        0.8125 0.93368 0.26983
        0.19132 0.93368 0.38577
        0.5 0.0 0.57733
        0.3125 0.75 0.81015
        0.19132 0.56632 0.38577
        0.8125 0.56632 0.26983
        0.0 0.25 0.46034
        0.0 0.0 0.07733
        0.75 0.25 0.67267
        0.80868 0.75 0.64719
        0.125 0.06632 0.82799
        0.125 0.43368 0.82799
        0.0 0.875 0.88708
        0.5 0.875 0.61293
        0.625 0.75 0.21034
        """

        self.coord= "relative"

        self.cages="""
        14 0.0 0.98471 0.32932
        14 -0.23471 -0.25 -0.07932
        12 -0.25 -0.25 -0.25
        14 0.5 -0.75 0.625
        14 0.0 0.48471 -0.32932
        12 0.5 1.5 1.0
        14 0.0 -0.98471 -0.32932
        12 0.5 0.75 0.08194
        12 0.25 -0.25 0.75
        15 0.5 0.25 0.52569
        14 1.0 1.51529 1.32932
        15 0.5 1.25 -0.27569
        12 0.0 1.0 0.5
        14 0.5 0.51529 0.17068
        14 0.0 -0.75 -0.125
        12 -0.25 0.25 -0.75
        15 0.0 -0.25 1.02569
        12 0.0 0.25 0.16806
        12 1.0 1.5 1.5
        12 0.0 -0.25 -0.16806
        14 -0.5 0.75 0.375
        14 0.5 0.98471 0.17068
        14 0.73471 0.75 0.57932
        12 0.5 1.0 1.0
        15 0.5 0.75 1.27569
        12 0.5 0.25 0.33194
        12 0.0 0.75 0.41806
        14 0.23471 -0.25 -0.07932
        14 0.23471 0.25 0.07932
        15 0.0 1.25 0.77569
        12 0.5 1.25 0.91806
        14 0.26529 0.25 0.42068
        15 0.0 0.25 -0.02569
        12 0.5 -0.25 0.66806
        14 0.76529 0.25 1.07932
        15 -0.5 -0.25 -0.52569
        15 0.0 -1.25 -0.77569
        14 0.5 1.01529 0.82932
        14 0.26529 -0.25 0.57932
        14 0.5 1.48471 0.82932
        14 -0.26529 0.25 -0.57932
        12 0.25 0.25 0.25
        12 0.0 -0.75 -0.41806
        14 0.0 0.75 0.125
        """

        self.bondlen = 3


        self.cell = """
        13.020869213334082 13.020869213334082 70.73623470390184
        """

        self.density = 0.6280734030039621



        self.cell = cellvectors(a=13.020869213334082,
                           b=13.020869213334082,
                           c=70.73623470390184)
