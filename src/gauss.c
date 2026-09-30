#include <math.h>
double genrand(void);

double gauss(void)
{
	static int have_extra=0;
	static double extra;
	double factor, hypot, x1, x2;

	if (have_extra == 0) {
		for (;;) {
			x1 = 2.0*genrand() - 1.0;
			x2 = 2.0*genrand() - 1.0;
			hypot = x1*x1 + x2*x2;
			if (hypot < 1.0 && hypot > 0.0) break;
		}

		factor = sqrt(-2.0*log(hypot)/hypot);
		extra = x1*factor;
		have_extra = 1;
		return x2*factor;
	} else {
		have_extra = 0;
		return extra;
	}
}
