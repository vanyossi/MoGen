#ifndef _MGN_TYPES_
#define _MGN_TYPES_

typedef struct _mgn_pop mgn_pop;
typedef struct _mgn_pop_param_pointer mgn_pop_param;
typedef enum _mgn_pop_type popType;

typedef struct _mgn_mop mgnMop;

typedef struct mgn_limit mgnLimit;

#define checkFlag(data, flag) (data & flag) == flag
#define UNUSED(x) ((void)(x))

#endif // _MGN_TYPES_
