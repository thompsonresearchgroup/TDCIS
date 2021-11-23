# TDCIS
A time-dependent CI solver that can use either an orthogonal or nonorthogonal determinant basis expansion.

The truncation level of the orthogonal CI expansion is arbitrary (singles only by default). A CAS expansion is not currelntly implemented, but is very easy to drop in. Nonorthogonal CI expansions currently require a single matrix file for each determinant in the expansion, which are typically obtained from a Gaussian calculation. 

Several different types of electric field pulses can be used:
1. Rectangular pulse
2. Delta pulse
3. Oscillating rectangular envelope
4. Oscillating Gaussian envelope
5. Chirped Gaussian envelope

Different pulse magnitudes, pulse lengths, pulse widths, chirp parameters, oscillation frequencies, and pulse onset/maximum time can be set using options.

Different inital states can be prepared.

Time correlation functions can be obtained from a specified initial time. 

Direct matrix elements can be computed without the requirement to store 2ERIs on disk (requires modified version of Gaussian).

See documentation folder or use the --help option for command line arguments used to specify options.

For any questions, please contact Lee M. Thompson, University of Louisville, lee.thompson.1@louisville.edu
(c) 2021 by Lee M. Thompson distributed under terms of MIT license.
