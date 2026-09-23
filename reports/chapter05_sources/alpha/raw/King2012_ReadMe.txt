J/MNRAS/422/3370    Spatial variation in fine-structure constant (King+, 2012)
================================================================================
Spatial variation in the fine-structure constant - new results from VLT/UVES.
    King J.A., Webb J.K., Murphy M.T., Flambaum V.V., Carswell R.F.,
    Bainbridge M.B., Wilczynska M.R., Koch F.E.
   <Mon. Not. R. Astron. Soc., 422, 3370-3414 (2012)>
   =2012MNRAS.422.3370K
================================================================================
ADC_Keywords: QSOs ; Redshifts ; Spectroscopy
Keywords: methods: data analysis - quasars: absorption lines -
          cosmology: observations

Abstract:
    Quasar absorption lines provide a precise test of whether the
    fine-structure constant, {alpha}, is the same in different places and
    through cosmological time. We present a new analysis of a large sample
    of quasar absorption-line spectra obtained using the Ultraviolet and
    Visual Echelle Spectrograph (UVES) on the Very Large Telescope (VLT)
    in Chile. We apply the many-multiplet method to derive values of
    {DELTA}{alpha}/{alpha}==({alpha}z-{alpha}0)/{alpha}0 from 154
    absorbers, and combine these values with 141 values from previous
    observations at the Keck Observatory in Hawaii.

Description:
    We have used publicly available spectra from VLT/UVES, in Chile (at
    25{deg} south), to generate an MM sample similar in size to the Keck
    sample.

File Summary:
--------------------------------------------------------------------------------
 FileName   Lrecl  Records   Explanations
--------------------------------------------------------------------------------
ReadMe         80        .   This file
tablea1.dat    56      295   Table of all VLT and Keck da/a values
--------------------------------------------------------------------------------

Byte-by-byte Description of file: tablea1.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   1-  3  I3    ---     Seq       Sequential number
   5- 18  A14   ---     QSO       QSO name (JHMMSS+DDMMSS)
  20- 24  F5.3  ---     zem       Emission redshift
  26- 32  F7.5  ---     zabs      Absorption redshift
  34- 40  F7.3  10-5    da/a      Value of ({Delta}{alpha})/{alpha}
                                  = ({alpha}_z_{alpha}_0_)/{alpha}_0_
                                  ({alpha}_z_ at the zabs redshift)
  42- 47  F6.3  10-5  e_da/a      Error on da/a
  49- 52  A4    ---     Sample    Sample (Keck or VLT)
      54  I1    ---     Flag      [1/3] sig_rand_flag value (1)
      56  I1    ---     Out       [0/1] indicates whether the system was
                                   identified as an outlier and removed from
                                   consideration (2)
--------------------------------------------------------------------------------
Note (1): The value of sigma_rand depends on the model under consideration.
  Some 11 models are considered in the paper, each with its own sigma_rand
  value. The values of sigma_rand used are given in Table 2. For a particular
  model, the sigma_rand values given in Table 2 should be added in quadrature
  with the "err" values according to the "sig_rand_flag" column:
  * sig_rand(VLT)     => sig_rand_flag = 3
  * sig_rand(Keck LC) => sig_rand_flag = 1
  * sig_rand(Keck HC) => sig_rand_flag = 2
  where LC and HC refer to "low contrast" and "high contrast"
  subsamples of the Keck sample. The particular case of a weighted
  mean requires:
  * sig_rand(VLT)     = 0.905
  * sig_rand(Keck LC) = 0
  * sig_rand(Keck HC) = 1.743
  [these values are in units of 10^-5^]
  There are 27 Keck "high contrast" (sig_rand_flag=2) absorbers.
  All the other Keck absorbers have sig_rand_flag=1, whilst
  all the VLT absorbers have sig_rand_flag=3.
Note (2): In the analysis in the paper, two outliers are identified through
  the Least Trimmed Squares analysis. These outliers were removed from
  consideration before any further statistical analysis. These can be
  identified by the presence of a "1" in the "Outlier" column (all
  other absorbers have "0" in the Outlier column).
--------------------------------------------------------------------------------

Table 2: sig_rand values used in each model.
--------------------------------------------------------------------------------
  I   Sample + model         sig_rand(VLT) sig_rand(Keck LC) sig_rand(Keck HC)
                               (10^-5^)      (10^-5^)         (10^-5^)
--------------------------------------------------------------------------------
  1   Keck04-dipole              N/A           0                1.630
  2   #1 with no monopole        N/A           0                1.668
  3   VLT-weighted mean          0.905         N/A              N/A
  4   Combined weighted mean     0.905         0                1.743
  5   VLT-dipole                 0.905         N/A              N/A
  6   #5 with no monopole        0.882         N/A              N/A
  7   Combined dipole            0.905         0                1.630
  8   #7 with no monopole        0.882         0                1.668
  9   Combined r-dipole          0.858         0                1.630
 10   #9 with no monopole        0.858         0                1.630
 11   z{beta} dipole             0.812         0                1.592
--------------------------------------------------------------------------------

History:
    From electronic version of the journal

================================================================================
(End)                                      Patricia Vannier [CDS]    24-Jan-2013
