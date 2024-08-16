#include <iostream>
#include "namespaces/settings.h"
#include "t2fit.h"
#include "tempo2.h"

void initialise_pulsar(int only_prefit)
{
    int num_pulsars = 1;

    std::cout << "form bats" << std::endl;

    formBatsAll(globals::pulsar, num_pulsars); /* Form Barycentric arrival times */
    std::cout << "form residuals" << std::endl;

    formResiduals(globals::pulsar, num_pulsars, 1); /* Form residuals */

    for (int it = 0; it < globals::num_tempo2_its; it++) {
        if (it > 0) /* Copy post-fit values to pre-fit values */
        {
            for (int i = 0; i < MAX_PARAMS; i++) {
                for (int k = 0; k < globals::pulsar->param[i].aSize; k++) {
                    globals::pulsar->param[i].prefit[k] = globals::pulsar->param[i].val[k];
                    globals::pulsar->param[i].prefitErr[k] = globals::pulsar->param[i].err[k];
                }
            }
        }

        for (int iteration = 0; iteration < 2; iteration++) /* Do pre- and post- fit analysis */
        {

            formBatsAll(globals::pulsar, num_pulsars);      /* Form Barycentric arrival times */
            formResiduals(globals::pulsar, num_pulsars, 1); /* Form residuals */

            if (iteration == 0) /* Only fit to pre-fit residuals */
            {
                t2Fit(globals::pulsar, num_pulsars, covarFuncFile);
            }

            if (iteration == 1 || only_prefit == 1) {
                /* Output results to the screen */
                textOutput(globals::pulsar, num_pulsars, 0, 0, 0, 0, 0);
            }
            globals::pulsar->noWarnings = 2;
            if (only_prefit == 1)
                iteration = 2;
        }
    }
}

void extra_delays(pulsar* psr, int num_pulsars);

void fastephemeris_routines(pulsar* psr, int num_pulsars)
{
    vectorPulsar(psr, num_pulsars);     /* 1. Form a vector pointing at the pulsar */
    readEphemeris(psr, num_pulsars, 0); /* 2. Read the ephemeris normally out */
    get_obsCoord(
        psr,
        num_pulsars); /* 3. Get Coordinate of observatory relative to Earth's centre normally out*/
    tt2tb(psr, num_pulsars);            /* Observatory/time-dependent part of TT-TB normally out*/
    readEphemeris(psr, num_pulsars, 0); /* Re-evaluate ephemeris with correct TB */
}

void fastformBatsAll(pulsar* psr, int num_pulsars)
{

    struct timeval tval_before, tval_after, tval_resultone;

    fastephemeris_routines(psr, num_pulsars); /* Ephemeris routines ... */

    extra_delays(psr, num_pulsars); /* Other time delays ... */

    formBats(psr, num_pulsars); /* Form Barycentric arrival times */

    secularMotion(psr, num_pulsars);
}