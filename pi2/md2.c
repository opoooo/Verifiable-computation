#include "md2.h"

int md2_init()
{
	// Set up the context
	if (core_init() != STS_OK)
	{
		core_clean();
		return 1;
	}
	if (pc_param_set_any() != STS_OK)
	{
		THROW(ERR_NO_CURVE);
		core_clean();
		return 1;
	}
	return 0;
}

int md2_close()
{
	core_clean();

	return 0;
}
