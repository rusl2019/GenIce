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
        183 8
        28 152
        92 195
        107 48
        18 197
        9 63
        133 150
        167 194
        69 169
        52 25
        57 212
        61 26
        216 82
        209 97
        171 176
        20 18
        59 199
        140 119
        113 40
        47 103
        226 73
        202 209
        27 131
        42 4
        66 82
        122 124
        76 104
        76 105
        130 209
        14 86
        97 54
        104 60
        77 106
        191 199
        15 92
        204 31
        165 67
        126 223
        128 213
        44 217
        157 206
        225 41
        168 103
        110 179
        111 177
        1 180
        26 210
        113 58
        28 222
        19 52
        135 217
        42 141
        77 62
        224 170
        32 0
        87 189
        204 55
        12 21
        102 150
        22 13
        146 35
        79 98
        13 23
        112 199
        39 147
        32 226
        21 144
        44 115
        154 143
        22 165
        58 200
        67 48
        125 227
        218 169
        140 11
        129 100
        178 48
        164 44
        212 109
        86 25
        180 115
        151 202
        88 17
        77 219
        45 205
        146 126
        19 190
        200 91
        93 2
        192 89
        160 138
        30 119
        208 53
        190 211
        157 12
        107 81
        121 174
        6 51
        225 222
        69 135
        39 175
        15 158
        131 132
        100 193
        149 50
        214 161
        215 203
        101 71
        121 134
        214 169
        90 112
        194 89
        81 119
        198 83
        193 10
        38 212
        84 106
        114 220
        129 208
        177 43
        107 173
        162 223
        166 98
        1 55
        134 141
        21 226
        155 126
        188 40
        172 196
        57 75
        20 181
        164 46
        3 27
        89 116
        120 41
        142 11
        129 9
        7 167
        18 142
        12 201
        208 125
        211 207
        101 37
        28 97
        5 124
        197 88
        80 126
        106 181
        183 30
        88 112
        90 196
        185 131
        19 135
        38 99
        43 98
        161 156
        185 174
        43 65
        38 220
        124 84
        117 179
        63 51
        37 62
        74 104
        121 148
        154 197
        37 173
        139 74
        110 166
        150 221
        27 168
        106 26
        32 11
        78 25
        154 101
        213 211
        207 16
        123 165
        157 49
        120 133
        159 55
        191 71
        158 56
        15 198
        130 221
        66 94
        23 180
        24 182
        33 88
        9 61
        104 51
        114 64
        100 130
        73 145
        201 172
        92 123
        208 163
        183 215
        73 168
        59 145
        0 145
        118 199
        50 91
        62 124
        207 159
        118 2
        143 26
        184 188
        184 187
        181 60
        218 156
        107 193
        114 117
        72 155
        131 47
        153 216
        142 63
        52 212
        181 136
        203 116
        195 127
        122 92
        86 224
        178 163
        192 172
        114 196
        41 97
        171 16
        134 144
        216 4
        75 179
        20 59
        113 156
        5 28
        223 69
        30 193
        170 209
        53 130
        217 85
        109 156
        94 87
        165 173
        46 162
        20 71
        117 36
        211 69
        197 210
        29 11
        172 167
        151 53
        21 45
        55 152
        127 222
        80 169
        149 205
        122 227
        29 74
        36 47
        5 56
        15 159
        162 217
        34 99
        10 83
        150 67
        100 54
        76 175
        96 123
        44 80
        140 206
        110 218
        190 86
        154 7
        8 10
        164 155
        31 176
        57 218
        40 144
        38 141
        187 65
        74 136
        61 56
        14 82
        12 34
        178 54
        3 34
        188 99
        129 6
        13 83
        153 72
        105 17
        207 195
        161 42
        60 210
        93 132
        128 195
        160 189
        125 61
        76 206
        7 182
        31 202
        121 65
        201 147
        75 220
        147 206
        202 152
        66 25
        57 135
        103 2
        168 187
        128 180
        167 186
        214 115
        191 89
        186 147
        153 176
        83 173
        151 227
        137 56
        75 85
        213 214
        101 194
        177 185
        3 174
        118 50
        95 65
        1 127
        34 134
        148 189
        196 91
        70 54
        79 47
        120 216
        170 78
        161 82
        98 146
        23 159
        93 91
        158 152
        123 62
        204 46
        203 81
        138 184
        158 227
        5 70
        63 210
        111 85
        148 4
        68 35
        95 146
        19 115
        166 80
        224 221
        0 33
        171 155
        29 203
        182 191
        220 200
        140 51
        143 219
        111 200
        139 6
        64 205
        145 2
        33 132
        10 221
        138 94
        138 95
        160 174
        108 136
        84 108
        148 72
        90 36
        87 78
        192 50
        14 213
        96 219
        110 58
        186 119
        184 3
        205 40
        171 127
        45 192
        102 41
        32 157
        13 224
        84 125
        96 8
        93 177
        105 73
        139 163
        226 175
        185 68
        24 186
        6 30
        111 166
        64 99
        120 170
        22 128
        1 164
        215 48
        182 219
        133 22
        37 137
        66 109
        18 0
        36 118
        29 175
        139 215
        223 16
        194 81
        201 64
        70 9
        24 183
        149 58
        122 222
        49 27
        149 117
        137 198
        43 103
        189 46
        87 31
        153 94
        49 105
        77 71
        132 90
        70 108
        79 68
        163 108
        160 35
        95 72
        204 16
        179 79
        176 225
        178 102
        142 136
        143 137
        49 33
        45 39
        7 112
        133 14
        162 35
        113 141
        53 102
        68 85
        109 188
        39 116
        96 198
        187 144
        78 4
        23 190
        24 116
        151 225
        60 17
        59 17
        8 67
        52 42
        """

        self.waters="""
        0.9501 0.5 0.92024
        0.4501 0.5 0.07976
        0.02495 0.681 0.03988
        0.13126 0.19392 0.75655
        0.33852 0.5 0.65205
        0.62996 0.5 0.96762
        0.75 0.125 0.66898
        0.86875 0.19392 0.24345
        0.6748 0.0 0.40575
        0.71245 0.319 0.84394
        0.63121 0.19392 0.452
        0.86868 0.5 0.72017
        0.02846 0.375 0.6212
        0.52846 0.375 0.3788
        0.44219 0.875 0.47171
        0.56876 0.125 0.1568
        0.40981 0.0 0.05176
        0.90981 0.0 0.94824
        0.87004 0.5 0.96762
        0.35474 0.319 0.34456
        0.8547 0.681 0.0355
        0.02846 0.625 0.6212
        0.52846 0.625 0.3788
        0.4846 0.319 0.28578
        0.82521 0.0 0.40575
        0.36879 0.19392 0.548
        0.75 0.125 0.02279
        0.06876 0.125 0.84321
        0.5499 0.5 0.92024
        0.85474 0.681 0.65544
        0.75 0.19392 0.54447
        0.43125 0.125 0.84321
        0.94875 0.5 0.76139
        0.97505 0.319 0.96012
        0.10988 0.31892 0.65165
        0.25 0.125 0.97721
        0.039 0.0 0.16681
        0.70506 0.5 0.22068
        0.18281 0.31892 0.46803
        0.94219 0.875 0.5283
        0.13121 0.80608 0.548
        0.51541 0.681 0.71422
        0.29358 0.5 0.52677
        0.14531 0.681 0.96451
        0.3251 0.5 0.1571
        0.0 0.75 0.5
        0.3547 0.319 0.96451
        0.09019 0.0 0.05176
        0.68281 0.68108 0.53198
        0.00351 0.194 0.83132
        0.01541 0.681 0.28578
        0.81116 0.194 0.73371
        0.31719 0.31892 0.46803
        0.5779 0.0 0.72859
        0.63132 0.5 0.72017
        0.47505 0.319 0.03988
        0.64531 0.319 0.0355
        0.25 0.125 0.33102
        0.14526 0.681 0.34456
        0.9172 0.806 0.02987
        0.84731 0.0 0.9094
        0.69024 0.194 0.95491
        0.68286 0.68108 0.16203
        0.78755 0.319 0.84394
        0.05781 0.125 0.47171
        0.18286 0.68108 0.83797
        0.32521 0.0 0.59425
        0.63121 0.80608 0.452
        0.19024 0.194 0.04509
        0.34871 0.0 0.2232
        0.6749 0.5 0.8429
        0.81715 0.68108 0.16203
        0.31715 0.68108 0.83797
        0.00351 0.806 0.83132
        0.81116 0.806 0.73371
        0.18885 0.194 0.26629
        0.92211 0.0 0.72859
        0.75 0.80608 0.14828
        0.39012 0.31892 0.65165
        0.15269 0.0 0.0906
        0.28755 0.681 0.15606
        0.79358 0.5 0.47323
        0.36879 0.80608 0.548
        0.60988 0.31892 0.34835
        0.69024 0.806 0.95491
        0.21245 0.319 0.15606
        0.44219 0.125 0.47171
        0.36875 0.19392 0.75655
        0.9172 0.194 0.02987
        0.89012 0.68108 0.34835
        0.99649 0.194 0.16868
        0.05126 0.5 0.23861
        0.56876 0.875 0.1568
        0.0499 0.5 0.07976
        0.32381 0.0 0.71165
        0.25 0.80608 0.85172
        0.67619 0.0 0.28835
        0.55126 0.5 0.76139
        0.19024 0.806 0.04509
        0.13121 0.19392 0.548
        0.64526 0.319 0.65544
        0.79494 0.5 0.22068
        0.57186 0.806 0.65137
        0.08281 0.806 0.97014
        0.84871 0.0 0.7768
        0.961 0.0 0.83319
        0.75 0.875 0.02279
        0.70642 0.5 0.47323
        0.71245 0.681 0.84394
        0.25 0.0 0.53591
        0.18885 0.806 0.26629
        0.1749 0.5 0.1571
        0.93125 0.125 0.1568
        0.18281 0.68108 0.46803
        0.07186 0.194 0.34864
        0.36868 0.5 0.27983
        0.86879 0.80608 0.452
        0.0779 0.0 0.27142
        0.99649 0.806 0.16868
        0.81719 0.31892 0.53198
        0.47154 0.625 0.6212
        0.20506 0.5 0.77932
        0.58281 0.806 0.02987
        0.63126 0.80608 0.24345
        0.64531 0.681 0.0355
        0.65269 0.0 0.9094
        0.30976 0.806 0.04509
        0.47505 0.681 0.03988
        0.4846 0.681 0.28578
        0.68885 0.194 0.73371
        0.57186 0.194 0.65137
        0.08281 0.194 0.97014
        0.02495 0.319 0.03988
        0.5 0.75 0.5
        0.16148 0.5 0.65205
        0.31116 0.194 0.26629
        0.78755 0.681 0.84394
        0.68286 0.31892 0.16203
        0.25 0.0 0.7707
        0.75 0.875 0.66898
        0.85474 0.319 0.65544
        0.20642 0.5 0.52677
        0.8251 0.5 0.8429
        0.75 0.19392 0.14828
        0.10988 0.68108 0.65165
        0.97505 0.681 0.96012
        0.25 0.875 0.97721
        0.94219 0.125 0.5283
        0.29494 0.5 0.77932
        0.07186 0.806 0.34864
        0.55781 0.875 0.5283
        0.539 0.0 0.83319
        0.52495 0.319 0.96012
        0.36875 0.80608 0.75655
        0.81715 0.31892 0.16203
        0.3547 0.681 0.96451
        0.25 0.80608 0.45553
        0.9846 0.319 0.71422
        0.58281 0.194 0.02987
        0.50351 0.194 0.16868
        0.25 0.19392 0.85172
        0.31719 0.68108 0.46803
        0.30976 0.194 0.04509
        0.68885 0.806 0.73371
        0.37004 0.5 0.03238
        0.60988 0.68108 0.34835
        0.21245 0.681 0.15606
        0.89012 0.31892 0.34835
        0.06876 0.875 0.84321
        0.31116 0.806 0.26629
        0.47154 0.375 0.6212
        0.4172 0.806 0.97014
        0.97154 0.375 0.3788
        0.66148 0.5 0.34795
        0.18286 0.31892 0.83797
        0.92814 0.806 0.65137
        0.43125 0.875 0.84321
        0.12996 0.5 0.03238
        0.64526 0.681 0.65544
        0.1513 0.0 0.2232
        0.44875 0.5 0.23861
        0.80976 0.806 0.95491
        0.82381 0.0 0.28835
        0.75 0.0 0.4641
        0.17619 0.0 0.71165
        0.14531 0.319 0.96451
        0.86879 0.19392 0.452
        0.13126 0.80608 0.75655
        0.1748 0.0 0.59425
        0.31715 0.31892 0.83797
        0.42814 0.194 0.34864
        0.86875 0.80608 0.24345
        0.97154 0.625 0.3788
        0.68281 0.31892 0.53198
        0.83852 0.5 0.34795
        0.50351 0.806 0.16868
        0.01541 0.319 0.28578
        0.8547 0.319 0.0355
        0.63126 0.19392 0.24345
        0.93125 0.875 0.1568
        0.13132 0.5 0.27983
        0.0 0.25 0.5
        0.49649 0.194 0.83132
        0.81719 0.68108 0.53198
        0.4172 0.194 0.97014
        0.05781 0.875 0.47171
        0.92814 0.194 0.65137
        0.461 0.0 0.16681
        0.6513 0.0 0.7768
        0.51541 0.319 0.71422
        0.80976 0.194 0.95491
        0.42211 0.0 0.27142
        0.25 0.19392 0.45553
        0.42814 0.806 0.34864
        0.35474 0.681 0.34456
        0.75 0.80608 0.54447
        0.39012 0.68108 0.65165
        0.28755 0.319 0.15606
        0.25 0.875 0.33102
        0.75 0.0 0.2293
        0.14526 0.319 0.34456
        0.55781 0.125 0.5283
        0.52495 0.681 0.96012
        0.34731 0.0 0.0906
        0.5 0.25 0.5
        0.49649 0.806 0.83132
        0.9846 0.681 0.71422
        0.59019 0.0 0.94824
        """

        self.coord= "relative"

        self.cages="""
        15 -0.5436 0.0 -0.3482
        16 0.07025 -0.5 -0.16699
        14 -0.42433 -0.5 -0.4684
        12 -0.4002 -0.27599 -0.15952
        12 -0.66117 0.0 -0.11199
        12 0.25 0.0 0.15534
        12 -0.0998 -0.27599 -0.15952
        15 -0.75 0.5 -0.02252
        14 0.42433 -0.5 0.4684
        12 0.4002 0.27599 0.15952
        12 0.0 0.0 0.0
        14 -0.25 -0.27567 -0.65339
        12 0.0998 -0.27599 0.15952
        15 0.75 0.5 0.02252
        15 0.25 0.5 0.33189
        14 0.25 0.27567 0.65339
        12 -0.15558 0.0 -0.41842
        12 0.15558 0.0 0.41842
        12 -0.25 0.0 -0.15534
        12 0.4002 -0.27599 0.15952
        15 -0.0436 0.0 0.3482
        12 0.0998 0.27599 0.15952
        12 0.66117 0.0 0.11199
        12 -0.16117 0.0 0.11199
        14 -0.07567 0.5 -0.4684
        12 -0.34442 0.0 -0.41842
        12 -0.0998 0.27599 -0.15952
        16 0.57025 0.5 0.16699
        12 0.34442 0.0 0.41842
        16 -0.07025 -0.5 0.16699
        15 -0.25 0.5 -0.33189
        14 -0.25 0.27567 -0.65339
        12 -0.4002 0.27599 -0.15952
        12 0.5 0.0 0.0
        12 0.16117 0.0 -0.11199
        14 0.07567 0.5 0.4684
        16 -0.57025 0.5 -0.16699
        15 0.0436 0.0 -0.3482
        14 0.25 -0.27567 0.65339
        15 0.5436 0.0 0.3482
        """

        self.bondlen = 3


        self.cell = """
        41.49002717228765 14.6803698879515 24.33425022050341
        """

        self.density = 0.45979888522587886



        self.cell = cellvectors(a=41.49002717228765,
                           b=14.6803698879515,
                           c=24.33425022050341)
