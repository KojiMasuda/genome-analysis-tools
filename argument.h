#ifndef _argument_h_INCLUDED_
#define _argument_h_INCLUDED_

typedef enum {
	ARGUMENT_TYPE_NONE,
	ARGUMENT_TYPE_FUNCTION,
	ARGUMENT_TYPE_FLAG_ON,
	ARGUMENT_TYPE_FLAG_OFF,
	ARGUMENT_TYPE_INTEGER,
	ARGUMENT_TYPE_FLOAT,
	ARGUMENT_TYPE_STRING
} Argument_Type;

#define ARGUMENT_FLAG_ON  1
#define ARGUMENT_FLAG_OFF 0

#define VER_MAJOR 1
#define VER_MOD 0
#define VER_MINOR 0

typedef struct _Argument {
	const char *name;
	Argument_Type type;
	void *value;
} Argument;

int argument_read (int *argcp, char *argv[], const Argument args[]);

#endif
