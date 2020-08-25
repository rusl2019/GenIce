# coding: utf-8
"""
Data source: Dutour Sikirić, Mathieu, Olaf Delgado-Friedrichs, and Michel Deza. “Space Fullerenes: a Computer Search for New Frank-Kasper Structures” Acta Crystallographica Section A Foundations of Crystallography 66.Pt 5 (2010): 602–615.

Cage composition:
 (12,14,15,16) = (20,12,0,8,)
"""

import genice.lattices
from genice.cell import cellvectors

class Lattice(genice.lattices.Lattice):
    def __init__(self):
        self.pairs="""
        93 202
        149 80
        157 89
        216 72
        48 177
        70 208
        223 119
        5 189
        6 189
        109 188
        224 71
        180 147
        69 44
        185 155
        120 32
        16 34
        19 190
        145 8
        137 223
        42 69
        153 208
        41 16
        37 178
        105 34
        86 59
        40 158
        100 97
        38 180
        183 60
        32 76
        101 86
        112 197
        53 181
        156 96
        37 204
        224 195
        135 147
        142 148
        219 192
        76 81
        11 162
        20 10
        11 158
        207 137
        41 201
        221 90
        30 182
        222 213
        41 168
        27 66
        218 85
        114 98
        54 122
        170 67
        135 50
        14 115
        150 71
        118 144
        0 91
        6 151
        13 74
        23 194
        194 57
        5 44
        90 56
        169 91
        113 75
        94 151
        59 225
        106 66
        142 203
        7 22
        48 72
        10 196
        109 96
        126 40
        68 214
        180 87
        153 11
        205 208
        16 214
        179 84
        154 82
        162 50
        187 158
        16 178
        126 203
        30 106
        160 200
        189 85
        32 102
        44 115
        18 35
        105 159
        123 144
        40 209
        101 39
        133 96
        68 105
        54 53
        31 30
        1 2
        127 126
        132 98
        174 128
        54 193
        118 140
        18 25
        149 182
        85 184
        214 171
        217 108
        7 202
        28 29
        23 78
        74 121
        0 41
        17 112
        20 114
        37 82
        201 36
        160 134
        155 128
        119 111
        153 79
        126 208
        58 77
        227 188
        133 63
        156 222
        157 226
        179 198
        95 168
        210 49
        200 175
        27 122
        64 108
        224 12
        49 164
        61 103
        136 36
        13 163
        14 165
        61 97
        35 164
        6 129
        160 31
        131 7
        22 204
        143 99
        161 4
        96 141
        15 207
        15 216
        112 211
        139 67
        8 105
        70 196
        69 46
        132 111
        11 90
        117 192
        156 51
        95 152
        151 184
        100 169
        217 163
        24 144
        37 130
        38 132
        75 177
        104 59
        139 165
        167 5
        32 95
        45 188
        170 73
        227 21
        68 104
        65 148
        211 183
        26 113
        101 68
        12 191
        127 157
        20 111
        103 52
        42 87
        199 34
        201 152
        172 93
        191 92
        3 12
        114 165
        58 176
        35 125
        117 212
        46 212
        190 186
        124 45
        29 123
        100 154
        42 89
        220 63
        19 170
        26 161
        138 43
        1 123
        134 23
        190 169
        120 143
        167 129
        187 203
        140 1
        53 141
        183 202
        204 183
        18 201
        194 197
        104 174
        26 136
        22 13
        153 212
        6 216
        213 55
        198 128
        135 42
        34 25
        59 22
        118 25
        196 161
        176 79
        8 39
        48 99
        176 89
        10 205
        119 146
        138 172
        206 155
        139 123
        125 166
        65 213
        215 65
        226 106
        17 52
        135 200
        216 75
        193 221
        149 213
        131 92
        94 113
        15 26
        124 138
        108 60
        206 159
        78 122
        176 46
        53 121
        192 161
        17 57
        2 177
        164 110
        175 141
        133 56
        129 219
        131 224
        131 225
        113 62
        58 78
        181 80
        88 103
        215 50
        75 14
        110 39
        173 211
        127 79
        7 45
        27 163
        84 82
        38 146
        52 3
        108 66
        199 150
        136 102
        33 152
        225 128
        12 57
        97 3
        185 100
        86 178
        69 215
        184 98
        205 209
        206 61
        179 88
        46 148
        9 73
        146 209
        81 154
        171 103
        209 147
        45 74
        17 179
        185 150
        33 102
        47 91
        15 33
        225 150
        174 159
        9 29
        0 171
        55 50
        55 175
        0 120
        19 49
        118 8
        130 52
        87 142
        20 4
        220 66
        106 107
        185 186
        35 33
        54 21
        117 70
        227 181
        173 13
        196 207
        62 116
        172 181
        64 191
        182 77
        125 137
        47 145
        144 110
        55 90
        117 44
        134 138
        24 28
        115 98
        172 74
        83 146
        210 39
        109 121
        145 206
        223 4
        67 164
        166 67
        93 173
        219 187
        31 51
        76 101
        119 158
        198 43
        72 116
        116 73
        160 107
        23 182
        27 121
        219 223
        222 56
        88 154
        10 83
        88 155
        38 85
        65 180
        212 187
        47 140
        71 97
        129 207
        221 141
        193 63
        84 217
        125 28
        9 166
        194 64
        140 19
        83 218
        30 64
        124 92
        197 122
        202 227
        217 204
        86 199
        190 99
        195 43
        226 77
        76 171
        49 25
        2 143
        188 51
        78 107
        168 81
        222 147
        218 166
        218 165
        14 189
        83 137
        28 4
        192 94
        93 195
        145 186
        156 200
        148 77
        167 70
        94 115
        220 21
        178 81
        130 214
        132 215
        18 24
        102 110
        82 3
        168 169
        139 177
        60 21
        5 87
        210 186
        184 9
        111 162
        95 99
        134 80
        114 29
        2 62
        159 71
        149 51
        47 143
        57 43
        112 163
        104 173
        31 124
        167 203
        175 80
        120 36
        60 197
        116 151
        127 133
        157 142
        210 199
        226 63
        220 109
        24 136
        48 170
        61 91
        130 211
        1 73
        92 198
        195 174
        36 62
        152 72
        79 221
        205 162
        40 56
        107 89
        58 193
        84 191
        """

        self.waters="""
        0.25585 0.02457 0.1875
        0.93406 0.9058 0.1875
        0.06594 0.9058 0.3125
        0.12844 0.17939 0.0
        0.5625 0.77399 0.1875
        0.12156 0.67174 0.625
        0.1875 0.76634 0.8125
        0.6875 0.26634 0.8125
        0.74415 0.02457 0.3125
        0.87156 0.82061 0.0
        0.625 0.71732 0.5
        0.56594 0.5942 0.1875
        0.03217 0.23748 0.0
        0.53217 0.26252 0.5
        0.03217 0.76252 0.5
        0.37844 0.82826 0.625
        0.47127 0.06691 0.0
        0.1875 0.23366 0.3125
        0.52873 0.93309 0.0
        0.90161 0.96369 0.875
        0.6875 0.73366 0.3125
        0.46651 0.37138 0.0
        0.5625 0.22601 0.6875
        0.05906 0.39076 0.1875
        0.55906 0.89076 0.1875
        0.65161 0.98716 0.0
        0.37844 0.82826 0.375
        0.37288 0.35362 0.5
        0.62156 0.82826 0.125
        0.77689 0.82444 0.1875
        0.05906 0.39076 0.8125
        0.94094 0.39076 0.6875
        0.3484 0.98716 0.5
        0.44094 0.89076 0.6875
        0.62288 0.05407 0.0
        0.55906 0.89076 0.8125
        0.31807 0.92136 0.1875
        0.37844 0.17174 0.875
        0.87844 0.67174 0.875
        0.65161 0.01284 0.5
        0.56594 0.5942 0.8125
        0.3484 0.01284 0.0
        0.02873 0.56691 0.5
        0.9375 0.27399 0.3125
        0.12156 0.67174 0.375
        0.72312 0.32444 0.6875
        0.18193 0.57864 0.1875
        0.0 0.0 0.25
        0.06594 0.9058 0.6875
        0.74415 0.97543 0.8125
        0.81807 0.57864 0.3125
        0.81807 0.42136 0.8125
        0.22312 0.17557 0.1875
        0.56594 0.4058 0.3125
        0.43406 0.4058 0.1875
        0.75585 0.52457 0.1875
        0.5984 0.53631 0.875
        0.0625 0.27399 0.1875
        0.24415 0.47543 0.1875
        0.62156 0.17174 0.625
        0.37156 0.32061 0.0
        0.06594 0.0942 0.1875
        0.19307 0.8752 0.1875
        0.40161 0.46369 0.875
        0.12156 0.32826 0.875
        0.97127 0.56691 0.0
        0.30693 0.3752 0.6875
        0.80693 0.8752 0.6875
        0.55906 0.10924 0.3125
        0.05906 0.60924 0.3125
        0.37156 0.67939 0.5
        0.87288 0.14639 0.0
        0.19307 0.8752 0.8125
        0.96651 0.87138 0.0
        0.62844 0.32061 0.5
        0.12844 0.82061 0.5
        0.37713 0.05407 0.5
        0.15161 0.48716 0.0
        0.18193 0.42136 0.3125
        0.40161 0.53631 0.375
        0.81807 0.42136 0.1875
        0.31807 0.07864 0.6875
        0.22312 0.17557 0.8125
        0.6875 0.73366 0.6875
        0.1875 0.23366 0.6875
        0.9375 0.72601 0.8125
        0.55906 0.10924 0.6875
        0.05906 0.60924 0.6875
        0.12713 0.14639 0.5
        0.15161 0.51284 0.5
        0.5984 0.53631 0.125
        0.0984 0.03631 0.125
        0.9375 0.27399 0.6875
        0.6875 0.26634 0.1875
        0.1875 0.76634 0.1875
        0.25585 0.97543 0.6875
        0.5984 0.46369 0.625
        0.03349 0.12862 0.0
        0.9375 0.72601 0.1875
        0.0984 0.96369 0.625
        0.06594 0.0942 0.8125
        0.52873 0.06691 0.5
        0.47127 0.93309 0.5
        0.19307 0.1248 0.3125
        0.62156 0.17174 0.375
        0.68193 0.07864 0.1875
        0.18193 0.42136 0.6875
        0.12288 0.44593 0.5
        0.27689 0.32444 0.8125
        0.56594 0.4058 0.6875
        0.62288 0.94593 0.5
        0.72312 0.67557 0.1875
        0.3125 0.26634 0.3125
        0.22312 0.82444 0.3125
        0.8125 0.76634 0.3125
        0.0625 0.72601 0.3125
        0.12713 0.85362 0.0
        0.27689 0.67557 0.3125
        0.74415 0.97543 0.1875
        0.62844 0.67939 0.0
        0.25585 0.97543 0.3125
        0.53349 0.37138 0.5
        0.30693 0.3752 0.3125
        0.80693 0.8752 0.3125
        0.87844 0.32826 0.625
        0.62156 0.82826 0.875
        0.43406 0.5942 0.6875
        0.40161 0.53631 0.625
        0.87156 0.17939 0.5
        0.3125 0.73366 0.8125
        0.37844 0.17174 0.125
        0.8125 0.23366 0.8125
        0.87844 0.67174 0.125
        0.5 0.5 0.75
        0.94094 0.39076 0.3125
        0.87713 0.55407 0.5
        0.44094 0.89076 0.3125
        0.5625 0.77399 0.8125
        0.87844 0.32826 0.375
        0.87288 0.85362 0.5
        0.90161 0.96369 0.125
        0.5984 0.46369 0.375
        0.18193 0.57864 0.8125
        0.0984 0.96369 0.375
        0.68193 0.92136 0.3125
        0.90161 0.03631 0.375
        0.72312 0.67557 0.8125
        0.81807 0.57864 0.6875
        0.12288 0.55407 0.0
        0.87713 0.44593 0.0
        0.80693 0.1248 0.8125
        0.125 0.78269 0.0
        0.31807 0.92136 0.8125
        0.43406 0.5942 0.3125
        0.19307 0.1248 0.6875
        0.96651 0.12862 0.5
        0.75585 0.47543 0.6875
        0.24415 0.52457 0.6875
        0.53349 0.62862 0.0
        0.80693 0.1248 0.1875
        0.97127 0.43309 0.5
        0.4375 0.77399 0.3125
        0.69307 0.6248 0.3125
        0.375 0.28269 0.5
        0.68193 0.92136 0.6875
        0.875 0.78269 0.5
        0.77689 0.82444 0.8125
        0.27689 0.67557 0.6875
        0.25585 0.02457 0.8125
        0.0984 0.03631 0.875
        0.93406 0.9058 0.8125
        0.31807 0.07864 0.3125
        0.72312 0.32444 0.3125
        0.5625 0.22601 0.3125
        0.77689 0.17557 0.3125
        0.75585 0.47543 0.3125
        0.24415 0.52457 0.3125
        0.03349 0.87138 0.5
        0.44094 0.10924 0.8125
        0.125 0.21732 0.5
        0.94094 0.60924 0.8125
        0.69307 0.3752 0.1875
        0.02873 0.43309 0.0
        0.46783 0.26252 0.0
        0.96783 0.76252 0.0
        0.93406 0.0942 0.6875
        0.90161 0.03631 0.625
        0.37288 0.64639 0.0
        0.69307 0.3752 0.8125
        0.0625 0.72601 0.6875
        0.0 0.0 0.75
        0.0625 0.27399 0.8125
        0.3125 0.73366 0.1875
        0.40161 0.46369 0.125
        0.12156 0.32826 0.125
        0.8125 0.23366 0.1875
        0.46783 0.73748 0.5
        0.27689 0.32444 0.1875
        0.96783 0.23748 0.5
        0.68193 0.07864 0.8125
        0.8484 0.48716 0.5
        0.37713 0.94593 0.0
        0.625 0.28269 0.0
        0.30693 0.6248 0.8125
        0.4375 0.22601 0.8125
        0.62713 0.64639 0.5
        0.93406 0.0942 0.3125
        0.4375 0.77399 0.6875
        0.46651 0.62862 0.5
        0.69307 0.6248 0.6875
        0.74415 0.02457 0.6875
        0.4375 0.22601 0.1875
        0.30693 0.6248 0.1875
        0.8484 0.51284 0.0
        0.44094 0.10924 0.1875
        0.94094 0.60924 0.1875
        0.22312 0.82444 0.6875
        0.3125 0.26634 0.6875
        0.8125 0.76634 0.6875
        0.375 0.71732 0.0
        0.43406 0.4058 0.8125
        0.5 0.5 0.25
        0.75585 0.52457 0.8125
        0.53217 0.73748 0.0
        0.875 0.21732 0.0
        0.77689 0.17557 0.6875
        0.24415 0.47543 0.8125
        0.62713 0.35362 0.0
        """

        self.coord= "relative"

        self.cages="""
        12 0.5 0.31537 0.25
        12 0.25 -0.25 0.5
        12 -0.5 -0.31537 -0.25
        14 0.62981 0.45306 0.0
        12 0.26377 0.09829 0.0
        16 0.12131 0.33066 0.5
        16 0.62131 0.16934 0.0
        14 0.0 0.5 -0.25
        16 0.12131 -0.33066 0.0
        14 0.12981 0.04694 0.5
        16 0.37869 0.83066 0.0
        14 -0.12981 -0.04694 0.5
        12 1.0 0.81537 -0.25
        12 -0.26377 0.09829 0.5
        12 -0.5 -0.31537 0.25
        12 0.23623 0.40171 0.0
        14 -0.12981 0.04694 0.0
        16 0.37869 0.16934 0.5
        14 0.5 0.0 0.25
        12 0.76377 0.40171 0.5
        12 0.25 0.25 0.0
        14 0.62981 0.54694 0.5
        12 0.23623 0.59829 0.5
        12 0.0 0.18463 -0.25
        16 0.62131 0.83066 0.5
        14 -0.5 0.0 -0.25
        12 0.5 0.31537 -0.25
        14 0.37019 0.54694 0.0
        14 0.37019 0.45306 0.5
        12 -0.25 0.25 0.5
        12 0.26377 -0.09829 0.5
        12 1.0 0.81537 0.25
        12 -0.26377 -0.09829 0.0
        14 0.12981 -0.04694 0.0
        12 0.0 0.18463 0.25
        16 -0.12131 -0.33066 0.5
        14 1.0 0.5 0.25
        12 -0.25 -0.25 0.0
        16 -0.12131 0.33066 0.0
        12 0.76377 0.59829 0.0
        """

        self.bondlen = 3


        self.cell = """
        18.875972990397965 44.75732628442065 13.339504157452888
        """

        self.density = 0.6047187466903077



        self.cell = cellvectors(a=18.875972990397965,
                           b=44.75732628442065,
                           c=13.339504157452888)
