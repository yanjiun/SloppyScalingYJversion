import os
import numpy
import Utils
reload(Utils)

"""
Useful information about the WindowScaling data sets
"""

# Type of simulation: Linear, NonLinear, KPZ etc
simulType = "KPZ"

# Rows to skip
# Attention:
# Data are load entirely, but the fitting is done
# skipping the first 'rows_to_skip' lines
rows_to_skip = 0
final_rows_to_skip =0

# Type of normalization for the function (None, NormBasic, NormIntegerSum, NormLog)
# None does not need to be a string, NormBasic and NormIntegerSum do...
normalization = None

# Decide to add corrections to scaling in the function module (True/False)
corrections_to_scaling = True

#how much to weight priors
priorCoeff = 1.0

# System Size (used to get the fileNames) and
# values of independent variables. (List of tuples if more than one)
#
#independentValues without window sizes...
#

Values = [[1024,[0.05, 0.1, 1,10]],\
          [4096,[0.001]]]

lambda_only_values = [[4096,[0.01,0.05,0.1,0.5,1]],\
                      [16384,[0.0016384,0.016384,0.16384]]]

smallest_kappa_lambda_only=[[16384,[0.0016384]]]

small_kappa_values = [[1024,[0.001]],[2048,[0.001]],[4096,[0.0001,0.001]]]

large_kappa_values = [[1024,[0.01,0.05,0.1,1]]]

smallest_kappa=[[8192,[0.08192]]]

linear_model_values =[[4096,[0.0004096,0.004096,0.04096,0.4096]]]

KPZ_values=[[4096,[0.0004096,0.004096,0.04096,0.4096,4.096,40.96]],\
            [8192,[0.0008192,0.008192,0.08192]],\
            [16384,[0.00016384,0.0016384,0.016384]]]

#KPZ_values_2= [[4096,[0.04096,0.4096,4.096,40.96]],\
#               [8192,[0.008192]],\
#               [16384,[0.016384]]]

#KPZ_values_single = [[4096,[0.4096]]]

#KPZ_values_single =[[4096,[40.96]]]

KPZ_values_single=[[16384,[0.0016384]]]

#non-spanning avs
KPZ_values_2=[[4096,[0.4096,4.096,40.96]],\
              [8192,[0.08192]],\
              [16384,[0.016384,0.0016384]]]

#Values=[\
#        [1024,[0.005]],\
#        [4096,[0.001,0.005]]\
#        ]

#Values=[\
#        [1024,[0.005]],\
#        [4096,[0.001,0.005,0.01,0.05,0.1,0.5,1]]\
#        ]

#Values = [\
#         [1024,[0.001,0.005,0.01,0.05,0.1,1,10]],\
#         [2048,[0.001,0.005]],\
#         [4096,[0.0001,0.001,0.005]]\
#         ]

#
#independentValues with window sizes
#
#List of curves with similar Ws
# Ws around 1.3
#independentValues=[(0.01, 1024, 128),(0.05, 1024, 64),(10,1024,8),(0.001,4096,512),(0.005,4096,256)] 
Ws1_3_list =  [\
            [1024,0.01,128],\
            [1024,0.05,64],\
            [1024,10,8],\
            [4096,0.001,512],\
            [4096,0.005,256]\
            ]
# Ws around 0.05
#independentValues = [(0.001, 1024, 16),(0.005, 1024, 8),(0.001,2048,16),(0.005, 2048,8),(0.0001,4096,64)]
Ws0_05_list =  [\
            [1024,0.001,16],\
            [1024,0.005,8],\
            [2048,0.001,16],\
            [2048,0.005,8],\
            [4096,0.0001,64]\
            ]
#Ws around 0.5
#independentValues = [(0.001,1024,128),(0.005,1024,64),(1,1024,8),(0.0001,4096,512)]
# Ws around 5
#independentValues = [(0.001, 1024, 1024),(0.005, 1024,512),(0.005, 4096, 1024),(0.001, 2048, 1024)]
# Ws LARGEST 30~170
#independentValues = [(1, 1024, 512),(10,1024,256),(1,1024,1024),(10,1024,512),(10,1024,1024)]

# for testing, large Ws and large k's
#independentValues = [(1,1024,256),(1,1024,512),(1,1024,1024),(10,1024,128),(10,1024,256),(10,1024,512),(10,1024,1024)]

# for testing, very small Ws...
#independentValues = [(0.001,1024,4),(0.001,2048,4),(0.001,4096,4)]

#independentValues = [(0.1,1024,8),(0.1,1024,16),(0.1,1024,32),(0.1,1024,64),(0.1,1024,128),(0.1,1024,256),(0.1,1024,512),(0.1,1024,1024)]


#List of curves with different Ws
#independentValues = [
#    (0.001, 2048, 8),(0.001,2048,256),
#    (0.005,1024,16),(0.05,1024,16),
#    (0.005, 1024,512),
#    (0.01,1024,64),(0.1,1024,64),
#    (0.005,4096,512)]
Ws_list = [\
            [1024,0.005,[16,512]],\
            [1024,0.05,16],\
            [1024,[0.01,0.1],64],\
            [2048,0.001,[8,256]],\
            [4096,0.005,512]\
            ]

#temp list of independentValues for jointModulerun (for making colors and symbols)
#independentValues = [(0.001, 1024, 16),(0.005, 1024, 8),(0.001,2048,16),(0.005, 2048,8),(0.0001,4096,64),(0.001, 2048, 8),(0.005,1024,16),(0.05,1024,16),(0.01,1024,64),(0.001,2048,256),(0.1,1024,64),(0.005,4096,512),(0.01, 1024, 128),(0.05, 1024, 64),(10,1024,8),(0.001,4096,512),(0.005,4096,256),(0.001,1024,1),(0.005,1024,1),(0.01,1024,1),(0.05,1024,1),(0.1,1024,1)]

# Using A00(system size) for A(s)
A00_list =  [\
            [1024,[0.001,0.005,0.01,0.05,0.1,1,10],1024],\
            [2048,[0.001,0.005],2048],\
            [4096,[0.0001,0.001,0.005],4096]\
            ]

# all A00 for KPZ (but for data sets without spanning avs)
# remove data sets with only 3 points for now (all W=4)
A00_KPZ_list = [\
                [4096,[4.096],8],\
                [4096,[4.096],(32,4096)],\
                [4096,[0.04096,0.4096,40.96],(8,4096)],\
                [8192,[0.008192,0.08192],(8,8192)],\
                [16384,[0.0016384,0.016384],(8,16384)]\
               ]

A00_KPZ_selection = [\
                      [4096,[0.4096,4.096],64],\
                      [4096,[0.4096],1024],\
                      [8192,[0.08192],64],\
                      [8192,[0.08192],1024],\
                      [8192,[0.08192],8192],\
                      [16384,[0.016384,0.0016384],512],\
                      [16384,[0.016384],2048]\
                    ]


A00_KPZ_no_large_kappa = [\
                      [8192,[0.08192],64],\
                      [8192,[0.08192],1024],\
                      [8192,[0.08192],8192],\
                      [16384,[0.016384,0.0016384],512],\
                      [16384,[0.016384],2048]\
                    ]



A00_KPZ_temp = [\
                 [4096,[4.096],64],\
                 [4096,[0.4096],1024],\
                 [8192,[0.08192],8192],\
                 [16384,[0.016384],2048],\
                 [16384,[0.0016384],512]]

A00_KPZ_debug=[[8192,[0.08192],1024]]

# Using A10's and A01's (sorted according to win/corr)
A10_list = [\
            #[1024,[0.001,0.005,0.01,0.005,0.01,0.05,0.1,10],(4,1024)],\
            [2048,[0.001,0.005],(4,2048)],\
            [4096,[0.0001,0.001,0.005],(4,4096)]\
            ]

# complete list for A10's and A01's KPZ
A10_KPZ_list =[\
               [4096,[0.04096,0.4096,4.096,40.96],(4,4096)],\
               [8192,[0.008192,0.08192],(4,8192)],\
               [16384,[0.0016384,0.016384],(4,16384)]\
               ]

A10_KPZ_selection=[\
                   [4096,[0.4096],64],\
                   [4096,[0.4096],256],\
                   [4096,[0.4096,4.096],2048],\
                   [8192,[0.08192],8],\
                   [8192,[0.08192],128],\
                   [8192,[0.08192],512],\
                   [16384,[0.016384,0.0016384],16],\
                   [16384,[0.016384],1024],\
                   [16384,[0.016384,0.0016384],8192]\
                  ]

A10_KPZ_no_large_kappa=[\
                   [8192,[0.08192],8],\
                   [8192,[0.08192],128],\
                   [8192,[0.08192],512],\
                   [16384,[0.016384,0.0016384],16],\
                   [16384,[0.016384],1024],\
                   [16384,[0.016384,0.0016384],8192]\
                  ]


A10_KPZ_debug =[[4096,[0.4096],256]]

#The COMPLETE relevant list of (L,k,W) for A11
A11_list = [\
            [1024,0.001, (1,512)],\
            [1024,[0.005,0.01],(1,256)],\
            [1024,[0.05,0.1],(1,128)],\
            [1024,1,(1,32)],\
            [1024,10,(1,16)],\
            [2048,0.001,(2,512)],\
            [2048,0.005,(2,256)],\
            [4096,0.0001, (2,2048)],\
            [4096,[0.001,0.005],(2,512)]\
            ]

#remove [(4096,4.096,256),(4096,0.0096,2048)] only 3 and 2 data points after binning
# The complete list for (L, k, W) for KPZ (A11)
A11_KPZ_list = [\
                [4096,[0.004096],(1,2048)],\
                [4096,[0.04096],(1,1024)],\
                [4096,[0.4096],(1,512)],\
                [4096,[4.096],(1,128)],\
                [4096,[40.96],(1,64)],\
                [8192,[0.008192,0.08192],(1,1024)],\
                [16384,[0.0016384],(1,4096)],\
                [16384,[0.016384],(1,4096)]\
                ]

A11_KPZ_Ws=[\
             [4096,[0.4096],(1,64)],\
             [4096,[40.96],(1,8)],\
             [4096,[4.096],(1,64)],\
             [16384,[0.0016384,0.016384],(1,1024)]\
             ]


# remove [16394,[0.0016384],8192] set for now
# replace with [16384,[0.016384],1024]
# remove[16384,[0.0016384],1]
A11_KPZ_selection =[\
                    [4096,[4.096],1],\
                    [16384,[0.0016384],1],\
                    [16384,[0.0016384],8],\
                    [16384,[0.0016384],64],\
                    [16384,[0.0016384],256],\
                    [16384,[0.0016384],1024],\
                    [16384,[0.016384],1024],\
                    [16384,[0.0016384],8192],\
                    [4096,[0.4096],512],\
                    [4096,[40.96],64]\
                    ]

A11_KPZ_no_large_kappa =[\
                    [16384,[0.0016384],1],\
                    [16384,[0.0016384],8],\
                    [16384,[0.0016384],64],\
                    [16384,[0.0016384],256],\
                    [16384,[0.0016384],1024],\
                    [16384,[0.016384],1024],\
                    [16384,[0.0016384],8192]\
                    ]


A11_W_1_KPZ =[\
                [4096,[0.0004096],1],\
                [4096,[0.004096,0.04096],1],\
                [4096,[0.4096],1],\
                [4096,[4.096],1],\
                [4096,[40.96],1],\
                [8192,[0.0008192,0.008192,0.08192],1],\
                [16384,[0.00016384],1],\
                [16384,[0.0016384],1],\
                [16384,[0.016384],1]\
                ]

A11_W_1 = [\
           [1024, [0.001,0.005,0.01,0.05,0.1,1,10], 1]
          ]
A11_W_2 = [\
           [1024, [0.001,0.005,0.01,0.05,0.1,1,10], 2]
           ]
A11_W = [\
           [1024, [0.001,0.005,0.01,0.05,0.1,1,10], 1]
           ]


#The COMPLETE list of (L, k, W)
Complete_list = [\
            [1024, [0.001,0.005,0.01,0.05,0.1,0.5,1,10],(1,1024)],\
            [2048, [0.001,005],(1,2048)],\
            [4096, [0.0001,0.001,0.005],(1,4096)]\
            ]

sortedValues = True
independentNames, independentValues = Utils.get_independent(A11_KPZ_selection, sorting = sortedValues)
#independentNames, independentValues = Utils.get_independent(A11_KPZ_list, sorting=sortedValues)
#independentNames, independentValues = Utils.get_independent(KPZ_values_2, sorting = sortedValues)
#indepdendentNamesA10,independentValuesA10 = Utils.get_independent(A10_KPZ_list,sorting=sortedValues)
independentNamesA10, independentValuesA10 = Utils.get_independent(A10_KPZ_selection, sorting=sortedValues)
#independentNamesA00, independentValuesA00 = Utils.get_independent(A00_KPZ_list, sorting=sortedValues)
independentNamesA00, independentValuesA00 = Utils.get_independent(A00_KPZ_selection,sorting = sortedValues)
independentNamesAsk,independentValuesAsk=Utils.get_independent(KPZ_values_2,sorting=sortedValues)

# Data to fit and models
#moduleNames = ['Ahk','Awk']
#moduleNames=['A11']
moduleNames = ['A11','A00','A10']
#moduleNames = ['A10', 'A11']
#moduleNames = ['Ahk','Ask','Awk']
#moduleNames = ['Awk','Ask']
#moduleNames = ['Awk','Ahk']
#moduleNames = ['Atant','Ahk']

# Directory where data to be fit is stored
if os.getlogin() == 'yj':
    #dataDirectory = "data/"
    #dataDirectory = "lambda_only_data/"
    dataDirectory = "KPZ_data_new/"
    #dataDirectory = "linear_model/"
    #dataDirectory = 'data_no_ones/'
elif os.getlogin() == 'gf':
    dataDirectory = "/home/meas/WinSim/NonLinear/data/"

Symbol, Color = Utils.MakeSymbolsAndColors(independentValues)
SymbolA10,ColorA10 = Utils.MakeSymbolsAndColors(independentValuesA10)
SymbolA00, ColorA00 = Utils.MakeSymbolsAndColors(independentValuesA00)
