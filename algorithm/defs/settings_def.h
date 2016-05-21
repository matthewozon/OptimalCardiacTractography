#ifndef _SETTINGS_DEF
#define _SETTINGS_DEF

//#define _OLD_GRAPH
#define HANDLE_VERTEX_MOVE //must be defined
#define VOXEL_CONSTRAINT 1

#define Pi 3.141592653589793238462643383279502884197169399375105820974944592L

//define how many edges can be moved at once
#define EDGES_PER_NODE 26
#define MAX_EDGE_MOVED_AT_ONCE 26

//sweep settings (for icm iteration)
#define RAND_SWEEP          1
#define DETERMINISTIC_SWEEP 0

//move settings
#define MOVE_T1    100
#define MOVE_T2    200
#define MOVE_T4    800
//#define MOVE_TALL 1600

#define MOVE_V    3200 //vertex move

//cooling settings
#define COOLING_EXP       0
#define COOLING_EXP_BLOCK 1

//energy settings
#define U1    10L
#define U15   50L ///scalar product normed over mean edge pairs + angle sigmoid
#define U18   80L ///scalar product normed over mean edge pairs + angle sigmoid (decreasing weight)

//#define CONCAVE_EXP 1

//SA settings
#define SA_EXP_COOLING        COOLING_EXP
#define SA_EXP_BLOCK_COOLING  COOLING_EXP_BLOCK

#define SA_U1  U1
#define SA_U13_2 U13_2
#define SA_U15 U15

#define SA_MOVE_1  MOVE_T1//100
#define SA_MOVE_2  MOVE_T2//200
#define SA_MOVE_12 MOVE_T1+MOVE_T2
#define SA_MOVE_4  MOVE_T4
#define SA_MOVE_V MOVE_V
#define SA_MOVE_V1 MOVE_V+MOVE_T1
#define SA_MOVE_V2 MOVE_V+MOVE_T2


//SC settings
#define SC_EXP_COOLING        COOLING_EXP
#define SC_EXP_BLOCK_COOLING  COOLING_EXP_BLOCK

#define SC_U1  U1
#define SC_U15 U15
#define SC_U18 U18

#define SC_MOVE_1_2 1300
#define SC_MOVE_1_2_4 159
#define SC_MOVE_V1 1250

#define SC_MOVE_V1V 1251 ///SC_MOVE_V1+ move only vertices at the end (same temperature estimation as SC_MOVE_V1)
#define SC_MOVE_V1PT 1252 ///SC_MOVE_V1+ quantity of vertex increase with iteration (temperature estimation: move init=SA_MOVE_1 end=SC_MOVE_V1)
//#define SC_MOVE_V1_START 1251
//#define SC_MOVE_V1_END 1252

//ICM settings
#define ICM_U1  U1
#define ICM_U15 U15

#define ICM_MOVE_1  MOVE_T1
#define ICM_MOVE_2  MOVE_T2

#define GRADUAL_ALPHA 1
#define FIXED_ALPHA   0

#define MARC_ICM    1
#define NO_MARC_ICM 0


//integration
#define MIDDLE_SUM       0
#define TRAPEZOIDAL_RULE 1
#define SIMPSON_RULE     2

//mask maker
#define NO_SAMPLING 0
#define NEAREST_NEIGHBOR 1
#define LINEAR 2
#define KERNEL_TRICK 4

//interpolation
//#define NEAREST_NEIGHBOR 1
#define LINEAR_INTERPOLATION LINEAR //2
#define GAUSSIAN_PHI 4
#define CAUCHY_PHI   8
#define INV_SQRT     16

#define VECTOR       10 ///will no longer exist
#define TENSOR       20
#define DWIS         40


//error handling
#define NULL_POINTER      1000
#define BAD_PARAM         2000
#define ENDED_NOT_SUCCESS 4000
#define _SUCCESS             0

#define DT_NONE                    0
#define DTI_UNKNOWN                0
#define DT_BINARY                  1
#define DT_UNSIGNED_CHAR          2
#define DT_SIGNED_SHORT           4
#define DT_SIGNED_INT              8
#define DT_FLOAT                   16
#define DT_COMPLEX                 32
#define DT_DOUBLE                  64
#define DT_RGB                     128
#define DT_ALL                     255


#define DICOM_ERROR 0
#define DICOM_BOOL 1
#define DICOM_CHAR 2
#define DICOM_UCHAR 4
#define DICOM_SHORT 8
#define DICOM_USHORT 16
#define DICOM_LONG 32
#define DICOM_ULONG 64
#define DICOM_FLOAT 128
#define DICOM_DOUBLE 256

//#define SMALL_NUM 1e-14
#define SMALL_NUM 1.7e-307
#define MAX_EDGES_PER_NODE 128

#endif
