%% OpenSeismoMatlab
% Software for strong ground motion data processing
%
% -
%% Help
% Help section of the OpenSeismoMatlab functions
%
% -
%%
% <help_baselineCorr.html baselineCorr.m>
%
% Baseline correction of acceleration time history
%
% -
%%
% <help_BLKIN.html BLKIN.m>
%
% Bilinear kinematic hysteretic model with elastic viscous damping
%
% -
%%
% <help_CDReSp.html CDReSp.m>
%
% Constant Ductility Response Spectra
%
% -
%%
% <help_CSReSp.html CSReSp.m>
%
% Constant Strength Response Spectra
%
% -
%%
% <help_DRHA.html DRHA.m>
%
% Dynamic Response History Analysis
%
% -
%%
% <help_FASp.html FASp.m>
%
% Fourier amplitude spectrum
%
% -
%%
% <help_ffft.html ffft.m>
%
% Fast Finite Fourier Transform
%
% -
%%
% <help_HalfStep.html HalfStep.m>
%
% Reproduce signal with half time step
%
% -
%%
% <help_IDA.html IDA.m>
%
% Incremental Dynamic Analysis
%
% -
%%
% <help_iffft.html iffft.m>
%
% Inverse Fast Finite Fourier Transform
%
% -
%%
% <help_LEReSp.html LEReSp.m>
%
% Linear Elastic Response Spectra
%
% -
%%
% <help_LIDA.html LIDA.m>
%
% Linear Implicit Dynamic Analysis
%
% -
%%
% <help_NLIDABLKIN.html NLIDABLKIN.m>
%
% Non Linear Implicit Dynamic Analysis of a BiLinear KINematic hardening
% hysteretic structure with elastic damping
% 
% -
%%
% <help_OpenSeismoMatlab.html OpenSeismoMatlab.m>
%
% Seismic parameters and processing of an acceleration time history
%
% -
%%
% <help_PulseDecomp.html PulseDecomp.m>
%
% Pulse decomposition of a velocity time history
%
% -
%%
% <help_RPSReSp.html RPSReSp.m>
%
% Rigid Plastic Sliding Response Spectra
%
% -
%% Examples
% Examples of various OpenSeismoMatlab applications
%
% -
%%
% <example_baselineCorr.html example_baselineCorr>
%
% Perform baseline correction
%
% -
%% 
% <example_CDReSp.html example_CDReSp>
%
% Calculate constant ductility response spectra of an artificially
% generated acceleration time history
%
% -
%%
% <example_Constant_ductility_response_spectra.html
% example_Constant_ductility_response_spectra>
%
% Calculate constant ductility response spectra of a suite of ground
% motions
%
% -
%%
% <example_Fourier_spectra.html example_Fourier_spectra>
%
% Calculate Fourier spectra of a suite of ground motions
%
% -
%%
% <example_general.html example_general>
%
% Test all functionalities of OpenSeismoMatlab in a single script
%
% -
%%
% <example_LEReSp.html example_LEReSp>
%
% Calculate linear elastic response spectra of an artificially
% generated acceleration time history
%
% -
%%
% <example_Linear_elastic_response_spectra.html
% example_Linear_elastic_response_spectra.m>
%
% Calculate the linear elastic response spectra of an earthquake suite
%
% -
%%
% <example_PulseDecomp.html PulseDecomp.m>
%
% Pulse decomposition of a velocity time history
%
% -
%%
% <example_RPSReSp.html RPSReSp.m>
%
% Calculate rigid plastic sliding response spectra
%
% -
%%
% <example_Spectra_comparison_1.html example_Spectra_comparison_1.m>
%
% Comparison of elastic and constant ductility response spectra for mu=1
%
% -
%%
% <example_Spectra_comparison_2.html example_Spectra_comparison_2.m>
%
% Comparison of constant ductility response spectra for mu=2
%
% -

%% Verification
% Validation of OpenSeismoMatlab results against the literature
%
%%
% <verification_CDRS.html verification_CDRS>
%
% Verification of a constant ductility response spectrum
%
% Reference: Chopra (2020)
%
% -
%%
% <verification_CSRS1.html verification_CSRS1>
%
% Verification of a constant strength response spectrum
%
% Reference: Tena-Colunga, A. (1999)
%
% -
%%
% <verification_CSRS2.html verification_CSRS2>
%
% Verification of a constant strength response spectrum
%
% Reference: Tena-Colunga, A. (1999)
%
% -
%%
% <verification_DRHA.html verification_DRHA>
%
% Calculate linear dynamic response of a MDOF shear building
%
% Reference: Chopra (2020)
%
% -
%%
% <verification_filter1.html verification_filter1>
%
% Verify the high pass Butterworth filter of OpenSeismoMatlab
%
% Reference: Boore (2005)
%
% -
%%
% <verification_filter2.html verification_filter2>
%
% Verify the low- and high- pass Butterworth filter of OpenSeismoMatlab
%
% Reference: Graizer (2012)
%
% -
%%
% <verification_Fourier_spectrum1.html verification_Fourier_spectrum1>
%
% Verify the Fourier amplitude spectrum of the San Fernando earthquake
% (1971) - Component N11E
%
% Reference: California Institute of Technology, EERL (1974)
%
% -
%%
% <verification_Fourier_spectrum2.html verification_Fourier_spectrum2>
%
% Verify the Fourier amplitude spectrum of the San Fernando earthquake
% (1971) - Component N79W
%
% Reference: California Institute of Technology, EERL (1974)
%
% -
%%
% <verification_IDA1.html verification_IDA1>
%
% Verify the incremental dynamic analysis for an acceleration time history
% with given duration D_5_75 and spectral acceleration
%
% Reference: Mashayekhi et al. (2020)
%
% -
%%
% <verification_IDA2.html verification_IDA2>
%
% Verify the incremental dynamic analysis for an earthquake suite
%
% Reference: Deng et al. (2017)
%
% -
%%
% <verification_IDA3.html verification_IDA3>
%
% Verify the incremental dynamic analysis for the Loma Prieta earthquake
% (1989) 
%
% Reference: Vamvatsikos & Cornell (2002)
%
% -
%%
% <verification_IDA4.html verification_IDA4>
%
% Verify the incremental dynamic analysis bias due to fitting of the
% capacity curve of a SDOF system with an elastoplastic bilinear fit
% according to FEMA-440
%
% Reference: De Luca et al. (2011)
%
% -
%%
% <verification_LIDA.html verification_LIDA>
%
% Calculate the dynamic response of a linear SDOF oscillator
%
% Reference: Chopra (2020)
%
% -
%%
% <verification_LIDA_NLIDABLKIN.html verification_LIDA_NLIDABLKIN>
%
% Compare the output of LIDA.m and NLIDABLKIN.m for a linear SDOF
% oscillator in terms of seismic input energy per unit mass.
%
% Reference: None
%
% -
%%
% <verification_NLIDABLKIN.html verification_NLIDABLKIN>
%
% Verify the energy time history of SDOF oscillator
%
% Reference: Chopra (2020)
%
% -
%%
% <verification_NLIDABLKIN2.html verification_NLIDABLKIN2>
%
% Verify the energy time histories of a nonlinear SDOF oscillator
%
% Reference: Uang, C. M., & Bertero, V. V. (1990)
%
% -
%%
% <verification_NLIDABLKIN3.html verification_NLIDABLKIN3>
%
% Verify the strength ductility relation of an elastoplastic system
%
% Reference: Mahin, S. A., & Lin, J. (1983)
%
% -
%%
% <verification_PulseDecomp1.html verification_PulseDecomp1>
%
% Verify the extraction of the velocity pulse of a strong motion velocity
% time history 
%
% Reference: Baker, J. W. (2007)
%
% -
%%
% <verification_PulseDecomp2.html verification_PulseDecomp2>
%
% Verify the extraction of the velocity pulses of a strong motion velocity
% time history and sum them to reconstruct the original ground motion
%
% Reference: Shahi, S. K., & Baker, J. W. (2014)
%
% -
%%
% <verification_SIH1952.html verification_SIH1952>
%
% Calculate the spectral intensity as defined by Housner (1952)
%
% Reference: Housner, G. W. (1952)
%
% -
%%
% <verification_SRS1.html verification_SRS1>
%
% Calculate the rigid plastic sliding response spectrum
%
% Reference: Garini, E., & Gazetas, G. (2016)
%
% -
%%
% <verification_SRS2.html verification_SRS2>
%
% Calculate the rigid plastic sliding response spectrum
%
% Reference: Paglietti, A., & Porcu, M. C. (2001)
%
% -

%% Copyright
%
% Copyright (c) 2018-2023 by George Papazafeiropoulos
%
% * Major, Infrastructure Engineer, Hellenic Air Force
% * Civil Engineer, M.Sc., Ph.D.
% * Email: gpapazafeiropoulos@yahoo.gr
%

