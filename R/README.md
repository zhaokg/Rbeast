# BEAST: Bayesian Change-Point Detection and Time-Series Decomposition

**BEAST** — the **Bayesian Estimator of Abrupt change, Seasonality, and Trend** — is a fast Bayesian model-averaging algorithm for decomposing time series and other one-dimensional sequential data into abrupt-change, trend, and seasonal/periodic components. BEAST is useful for:

- change-point detection: breakpoints, structural breaks, joinpoints, regime shifts, and anomalies;
- trend analysis and nonlinear trend detection;
- decomposition of trend and seasonal/periodic components (i.e., seasonal-Trend decomposition
);
- time-series segmentation;
- interrupted time-series analysis;
- outlier detection;
- curve fitting and smotthing
- gap-filing of 1D curves

The algorithm is described in [Zhao et al. (2019)](https://drive.google.com/file/d/1MFZ0FpK1NwTieVSAf5jicLgl85Lm48uh/view). Examples on the use of BEAST/Rbeast across remote sensing, ecology, hydrology, public health, finance, paleoclimate, and other fields are provided in [**<ins>selected publications using BEAST/Rbeast</ins>**](#selected-publications-using-beastrbeast).

---

## Contents

- [Quick installation](#quick-installation)
- [Quick examples](#quick-examples)
- [Installation](#installation)
  - [R](#r-id)
  - [Python](#python-id)
  - [MATLAB](#matlab-id)
  - [Octave](#octave-id)
- [Description of BEAST](#description-of-beast)
- [Notes on computation](#notes-on-computation)
- [Compilation from C source code](#compilation-from-c-source-code)
- [References](#references)
- [Reporting bugs and requesting help](#reporting-bugs-and-requesting-help)
- [Acknowledgements](#acknowledgements)
- [Selected publications using BEAST/Rbeast](#selected-publications-using-beastrbeast)

---

## Quick installation

BEAST is implemented in C/C++ as a package named `Rbeast` and is available through R, Python, MATLAB, and Octave interfaces.

| Interface | Installation command |
|---|---|
| Python | `pip install Rbeast` |
| R | `install.packages("Rbeast")` |
| MATLAB | `eval(webread('http://b.link/rbeast', weboptions('cert','')))` |
| Octave | `eval(webread('http://b.link/rbeast'))` |

---

## Quick examples

### Python

```python
import Rbeast as rb

nile, year = rb.load_example("nile")
out = rb.beast(nile, start=1871, season="none")
rb.print(out)
rb.plot(out)
```

### MATLAB / Octave

```matlab
load('Nile.mat')
out = beast(Nile, 'start', 1871, 'season', 'none');
printbeast(out)
plotbeast(out)
```

### R

```r
library(Rbeast)

data(Nile)
out <- beast(Nile, season = "none")
print(out)
plot(out)
```

---

# Installation

## R  <a name=r-id> </a>  [![CRAN version](https://www.r-pkg.org/badges/version/Rbeast?color=green)](https://cran.r-project.org/package=Rbeast)   [![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/Rbeast?color=green)](https://cran.r-project.org/package=Rbeast)  [![CRAN Task View](https://img.shields.io/static/v1?style=plastic&logo=r&label=Rbeast&message=In%20CRAN%20Task%20Views&color=brightgreen)](https://cran.r-project.org/package=Rbeast)

`Rbeast` is available from CRAN. It is also listed in several CRAN Task Views, including [Time Series Analysis](https://cran.r-project.org/web/views/TimeSeries.html#forecasting-and-univariate-modeling), [Bayesian inference](https://cran.r-project.org/web/views/Bayesian.html#time-series-models), and [Environmetrics](https://cran.r-project.org/web/views/Environmetrics.html#environmental-time-series). Install it from CRAN:

```r
install.packages("Rbeast")
```

 **Note:** CRAN also hosts another package named `beast`, which is unrelated to this project. The package described here is `Rbeast`. It is also unrelated to the evolutionary-analysis software BEAST, which stands for Bayesian Evolutionary Analysis by Sampling Trees.

### Run and test Rbeast in R

The main R functions are `beast()`, `beast.irreg()`, and `beast123()`.

```r
library(Rbeast)

data(Nile)                          # Annual streamflow of the Nile River
out <- beast(Nile, season = "none") # Trend-only data without seasonality
print(out)
plot(out)

?Rbeast                            # See package documentation
```

<table border="0"  style='border:none;'  bordercolor="#ffffff"  width=100%  >
<tr style='border:none;'  >   
   <td valign="center" style='border:none;'  > 
        <img  height="300" align="left" src="https://github.com/zhaokg/Rbeast/raw/master/R/Images/beach.png">
   </td>
   <td valign="center" style='border:none;'  > 
       <img  height="300" align="center"  src="https://github.com/zhaokg/Rbeast/raw/master/R/Images/Nile.png">
   </td>   
 </tr>
 </table>

The package also includes two small games:

```r
tetris()                         # if you dare to waste a few moments of your life 
minesweeper()                    # if you dare to waste a few more moments of your life 
```

---

## Python   <a name=python-id></a>   [![PyPI](https://img.shields.io/static/v1?style=plastic&logo=python&label=Python&message=pip%20install%20Rbeast&color=brightgreen)](https://pypi.org/project/Rbeast/)

`Rbeast` is available from PyPI at  https://pypi.org/project/Rbeast/. Install the binary wheel package using:

```bash
pip install Rbeast
```

Binary wheel files are available for Windows, macOS, and Linux for common Python versions and CPU architectures. If the installation above from a wheel fails, install from source:

```bash
pip install Rbeast --no-binary :all:
```

> **Note**: Building from source requires a C/C++ compiler, such as MinGW GCC on Windows, GNU GCC on Linux or Xcode/Clang on macOS. If needed, contact Kaiguang Zhao (zhao.1423@osu.edu) to help build the package for your OS platform and Python version.

### Run and test Rbeast in Python

`Nile` is the annual streamflow of the River Nile, starting in 1871. Because these are annual observations, the series has no seasonal component.

```python
import Rbeast as rb

nile, year = rb.load_example("nile")
out = rb.beast(nile, start=1871, season="none")
rb.print(out)
rb.plot(out)

out  # Show the output fields
```

The second example, `googletrend`, is a monthly time series of Google Search popularity for the word **beach** in the United States. It is regularly spaced and has a yearly periodic component. Since the time step is one month, there are 12 data points per year.

```python
import Rbeast as rb

beach, year = rb.load_example("googletrend")

# Equivalent ways to specify a monthly time step and a yearly period
out = rb.beast(beach, start= 2004.0, deltat=1/12, period = 1.0)       # the time unit is unknown or arbitrary
out = rb.beast(beach, start= 2004.0, deltat=1/12, period ='1.0 year') # the time unit is fractional year
out = rb.beast(beach, start= 2004.0, deltat='1 month', period =1.0)   # the time unit is fractional year

rb.print(out)
rb.plot(out)
```

---

## MATLAB    <a name=matlab-id> </a>  [![View Rbeast on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/72515-bayesian-changepoint-detection-time-series-decomposition)

Install the MATLAB version of BEAST automatically to a local folder of your choice:

```matlab
beastPath = 'C:\beast';                 % Specify a target folder
eval(webread('http://b.link/rbeast'))   % Install to beastPath
```

If `webread` gives a certificate error, try:

```matlab
beastPath = 'C:\beast';
eval(webread('http://b.link/rbeast', weboptions('cert','')))
```

**Notes**:  <br/>
 \- The above will automatically download the files in the `Rbeast\Matlab` folder at Github to the chosen local path <br/>
 \-  You need write permission for the folder specified by `beastPath`.<br/>
 \-  The variable name must be exactly `beastPath`.<br/>
 \- If `beastPath` is not specified, the installer uses a temporary folder by default.<br/>
 \- If automatic installation fails, manually download the MATLAB files from the [Rbeast GitHub repository](https://github.com/zhaokg/Rbeast).<br/>

The MATLAB distribution includes:

- a compiled MEX binary bary library from the C/C++ soure code, such as `Rbeast.mexw64` on Windows, `Rbeast.mexa64` on Linux, or `Rbeast.mexmaci64`/`Rbeast.mexmaca64` on macOS;
- MATLAB wrapper functions such as `beast.m` and `beast123.m`;
- example datasets such as `Nile.mat` and `co2.mat`.

Precompiled MEX binaries such as Rbeast.mexw64 and Rbeast.mexa64 are provided for  Windows, Linux, and MacOS. If the included MEX binary does not work on your machine, you can compile it from the C source files in the [`Rbeast/Source`](https://github.com/zhaokg/Rbeast/tree/master/Source) folder. If needed, we are happy to work with you to compile for your specific machine. Additional information on compilations from the C source is also given below.

### Run and test Rbeast in MATLAB

The MATLAB  API is similar to those of R. Below is a quick example:

```matlab
help beast
help beast123

load('Nile.mat')
out = beast(Nile, 'season', 'none', 'start', 1871);
printbeast(out)
plotbeast(out)
```

---

## Octave  <a name=octave-id> </a> 

The Octave interface is similar to the MATLAB interface. Currently, the precompiled Octave version is primarily supported on Windows. For Octave on Linux or macOS, please contact Kaiguang Zhao at <zhao.1423@osu.edu> for assistance.

```octave
eval(webread('http://b.link/rbeast'))
```

---

## Julia and IDL  

Wrappers for Julia and IDL are under development. Contributions are welcome. Interested developers may contact Kaiguang Zhao at <zhao.1423@osu.edu>.

---

## Description of BEAST

Interpreting time-series data is often affected by model choice. Different models can produce different, or even contradictory, estimates of trends, patterns, and mechanisms from the same data. BEAST addresses this limitation by moving away from a single-best-model strategy and using Bayesian model averaging across competing models.

BEAST is designed to detect abrupt changes, cyclic or seasonal variations, and nonlinear trends in time-series observations. It not only estimates when changes occur but also quantifies the probability that each detected change is real. It can detect piecewise linear trends as well as more flexible nonlinear trends.

BEAST is applicable to many types of real-valued time-series or sequential data, including data from remote sensing, finance, public health, economics, climate science, ecology, hydrology, bioinformatics, education, sociology, and other fields. Example applications include identifying ecological regime shifts, mapping forest disturbance and land degradation from satellite imagery, detecting market trends, identifying anomalies and extreme events in climate data, and studying system dynamics in biological data.

Details are provided in [Zhao et al. (2019)](https://drive.google.com/file/d/1MFZ0FpK1NwTieVSAf5jicLgl85Lm48uh/view). The paper is also available at <https://go.osu.edu/beast2019>.

---

## Notes on computation

BEAST is a Bayesian algorithm, but it is designed to be computationally efficient--possibly among the fastest implementations of Bayesian time-series analysis algorithms of the same nature. (However, compared to non-Bayesian methods, it may be still slower due to the MCMC sampling-based inference.) 

For applications involving a few to thousands of time series, computation is usually manageable and not a practical concern. For remote-sensing and geospatial applications involving millions or billions of time series, computation can become a major challenge on desktop computers. For large stacked time-series images, first test BEAST on a single time series or a small image chip. This helps determine whether BEAST is appropriate for the application and provides a practical estimate of total processing time.

For 3D stacked image time series, use **`beast123`** rather than `beast` or `beast.irreg`. The `beast123` function is designed for data cubes and supports internal parallel computing. We also welcome consultation with Kaiguang Zhao (zhao.1423@osu.edu) to give specific suggestions if you see some value of BEAST for your applications.

---
 

## Compilation from C source code

Most users do not need to compile BEAST manually. However, developers and advanced users can compile the code for specific machines or platforms. The C/C++ source files are available in the [`Rbeast/Source`](https://github.com/zhaokg/Rbeast/tree/master/Source) folder.

For more detailed platform-specific compilation commands, see the source-code compilation instructions in the repository.

---

## References

- Zhao, K., Wulder, M. A., Hu, T., Bright, R., Wu, Q., Qin, H., Li, Y., Toman, E., Mallick B., Zhang, X., & Brown, M. (2019). [Detecting change-point, trend, and seasonality in satellite time series data to track abrupt changes and nonlinear dynamics: A Bayesian ensemble algorithm.](https://drive.google.com/file/d/1MFZ0FpK1NwTieVSAf5jicLgl85Lm48uh/view) Remote Sensing of Environment, 232, 111181. (the BEAST paper) 

- Zhao, K., Valle, D., Popescu, S., Zhang, X. and Mallick, B., 2013. [Hyperspectral remote sensing of plant biochemistry using Bayesian model averaging with variable and band selection](https://drive.google.com/file/d/1WGOAbH_h5f8ptQvsJsZqwlykO2JevsgT/view?usp=sharing). Remote Sensing of Environment, 132, pp.102-119. (the mcmc sampler used for BEAST)

- Hu, T., Toman, E.M., Chen, G., Shao, G., Zhou, Y., Li, Y., Zhao, K. and Feng, Y., 2021. [Mapping fine-scale human disturbances in a working landscape with Landsat time series on Google Earth Engine](https://pages.charlotte.edu/gang-chen/wp-content/uploads/sites/184/2021/05/Hu_2021_BEAST-HF-s.pdf). ISPRS Journal of Photogrammetry and Remote Sensing, 176, pp.250-261. (an application paper)

---

## Reporting bugs and requesting help

BEAST is distributed as-is, without warranty. Although it has been widely used, bugs and platform-specific issues often occur.

To report a bug, please include:

- your operating system;
- the interface used: R, Python, MATLAB, or Octave;
- a minimal reproducible example, if relevant;
- the full error message.

You may report issues through the [GitHub Issues tracker](https://github.com/zhaokg/Rbeast/issues) or contact Kaiguang Zhao at <zhao.1423@osu.edu>.  We also welcome any inquiry on feature enhancements or specific usage advice.
 

---

## Acknowledgements

BEAST was developed by Yang Li, Tongxi Hu, Xuesong Zhang, and Kaiguang Zhao.  Development of BEAST were partly funded by Microsoft Azure for Research (CRM0518513), a USGS 104B grant, and a Harmful Algal Bloom Research Initiative grant from the Ohio Department of Higher Education. Xuesong Zhang’s contribution was supported by USDA-ARS.

---

 

## Selected publications using BEAST/Rbeast

BEAST was originally published for ecological and environmental applications, but it is a generic statistical tool for time-series and sequential data across many disciplines. BEAST is not a heuristic algorithm but a rigorous statistical model. To illustrate the breadth of its applications, below are some selected pulications that used BEAST for statistical data analysis: 

| Discipline | Publication Title |
| --- | --- |
| Remote Sensing| *Li, J., Li, Z., Wu, H., and You, N., 2022. [Trend, seasonality, and abrupt change detection method for land surface temperature time-series analysis: Evaluation and improvement](https://doi.org/10.1016/j.rse.2022.113222). Remote Sensing of Environment, 10.1016/j.rse.2022.113222*|
|Archaeogenetics|*Sikora, M., Canteri, E., Fernandez-Guerra, A., Oskolkov, N., Ågren, R., Hansson, L., Irving-Pease, E.K., Mühlemann, B., Holtsmark Nielsen, S., Scorrano, G. and Allentoft, M.E., 2025. The spatiotemporal distribution of human pathogens in ancient Eurasia. Nature, 643(8073), pp.1011-1019.*|
|Energy|*Jakhmola, A., Jewell, J., Vinichenko, V. et al.2026. Probabilistic projections of global wind and solar power growth based on historical national experience. Nature Energy. https://doi.org/10.1038/s41560-026-02021-w*|
|Sustainability|*Sun, X., Tian, L., Fang, H., Walling, D.E., Huang, L., Park, E., Li, D., Zheng, C. and Feng, L., 2025. Changes in global fluvial sediment concentrations and fluxes between 1985 and 2020. Nature Sustainability, 8(2), pp.142-151.*|
|Paleoclimatology|*Lu, F., Lu, H., Gu, Y., Lin, P., Lu, Z., Zhang, Q., Zhang, H., Yang, F., Dong, X., Yi, S. and Chen, D., 2025. Tipping point-induced abrupt shifts in East Asian hydroclimate since the Last Glacial Maximum. Nature Communications, 16(1), p.477.*|
|Biogeochemistry|Sun, X., Tian, L., Fang, H., Walling, D.E., Syvitski, J., Huang, L., Li, D., Zheng, C. and Feng, L., 2026. Mapping pan-Arctic riverine particulate organic carbon from space (1985 to 2022). Science Advances, 12(3), p.eady6314.|
|Global Ecology|*Guo, R., Wu, X., Wang, P., Chen, T., Chen, X., Cai, J., Wang, X., Zhang, Z., Meng, Z. and Liu, Y., 2026. Increased spread of global flash droughts threatens vegetation productivity resilience. Nature Communications.*|
|Paleoceanography|*Rahaman, W., Gutjahr, M. and Prabhat, P., 2025. Late Pliocene growth of the West Antarctic Ice Sheet to near-modern configuration. Nature Communications, 16(1), p.6705.*|
|Zoology|*Meireles, J.P., Hahn-Klimroth, M., Lackey, L.B., van Eeuwijk, N., Bertelsen, M.F., Dressen, S., Dierkes, P.W., Abraham, A.J. and Clauss, M., 2026. Aging populations threaten conservation goals of zoos. Proceedings of the National Academy of Sciences, 123(5), p.e2522274123.*|
|Medicial sciences|*Garcia-Aymerich, J., de Las Heras, M., Carsin, A.E., Accordini, S., Agustí, A., Bui, D., Dharmage, S.C., Dodd, J.W., Eze, I., Gehring, U. and Gislason, T., 2025. General population-based lung function trajectories over the life course: an accelerated cohort study. The Lancet Respiratory Medicine, 13(7), pp.611-622.*|
|Climate Change|*Leclercq, L., Oelsmann, J., Cazenave, A., Passaro, M., Jevrejeva, S., Connors, S., Legeais, J.F., Birol, F. and Abarca-del-Rio, R., 2026. Abrupt trend change in global mean sea level and its components in the early 2010s. Communications Earth & Environment.*|
|Communication|*Phillips, J.B., 2025. Exploring psyop-based conspiracy theories on social media. Information, Communication & Society, pp.1-20.*|
|Applied Math|*Koutrouli, E., Manousopoulos, P., Theal, J. and Tresso, L., 2025. Crypto asset markets vs. financial markets: Event identification, latest insights and analyses. AppliedMath, 5(2), p.36.*|
|Biosystems Engineering|*Mayrhuber, E., Maschat, K., Brunner, D., Winkler, S.M. and Oczak, M., 2026. Improved and interpretable accelerometer-based farrowing prediction. Biosystems Engineering, 263, p.104381.*|
|Finance|*Li, S., Mishra, T. and Yarovaya, L., 2026. Heterogeneous Market Efficiency in Cryptocurrency Markets: A Multi‐Frequency Memory‐Based Approach. International Journal of Finance & Economics.*|
|Economics|*Lakštutienė, A., Sutiene, K., Kabasinskas, A., Malakauskas, A. and Kopa, M., 2025. Sustaining in Uncertain Time: Investigating Pension Fund Performance during Market Stress. Engineering Economics, 36(1), pp.96-112.*|
|Reliability Engineering|*Zhou, Y., Liu, S., Kou, G. and Kang, F., 2025. Degradation variation pattern mining based on BEAST time series decomposition integrated functional principal component analysis. Reliability Engineering & System Safety, 259, p.110952.*|
|Biology|*Bancaud, A., Nakajima, T., Suehiro, J.I., Alric, B., Morfoisse, F., Cacheux, J. and Matsunaga, Y.T., 2025. Intraluminal pressure triggers a rapid and persistent reinforcement of endothelial barriers. Lab on a Chip, 25(8), pp.2061-2072.*|
|Construction Engineering|*Li, F., Xie, Z., Yu, X. and Shi, B., 2025. Estimation of long-term variation patterns in the modal properties of a skyscraper under environmental effects. Engineering Structures, 336, p.120451.*|
|Animal Ecology|*Hole, G.M., Büntgen, U., Wang, Y., DeVries, B., Rees, G. and Wheeler, H.C., 2026. Dendrochronology and remote sensing reveal beaver occupancy and colonization dynamics in an expanding Arctic population. Ecosphere, 17(3), p.e70557.*|
|Wildlife|*Pérez-Mellado, V. and Pérez-Cembranos, A., 2025. Effects of the introduction of an herbivore on an endangered lizard. European Journal of Wildlife Research, 71(3), p.47.*|
|Coastal & Estuary science |*Richey, A., Oberbauer, S.F., Castaneda-Moya, E., Troxler, T., Kominoski, J.S., Olivas, P. and Malone, S.L., 2025. Sea-level rise and freshwater management are reshaping coastal landscapes. Journal of Environmental Management, 387, p.125842.*|
|Civil Engineering|*Englezou, Y., Timotheou, S. and Panayiotou, C.G., 2025. Fault-adaptive traffic demand estimation using network flow dynamics. IEEE Transactions on Intelligent Transportation Systems.*|
|Computational Chemistry|*Faran, M. and Bisker, G., 2025. Coarse-graining self-assembly by the stochastic landscape method. Journal of Chemical Theory and Computation, 21(21), pp.10719-10734.*|
|Education|*Huang, X., Wu, H., Liu, X. and Lajoie, S.P., 2025, July. What makes teamwork work? A multimodal case study on emotions and diagnostic expertise in an intelligent tutoring system. In International Conference on Artificial Intelligence in Education (pp. 44-52). Cham: Springer Nature Switzerland.*|
|Agrometerology|*Ssembajwe, R., Mulinde, C., Ddumba, S.D., Kagezi, G.H., Opio, R., Kobusinge, J., Mugagga, F., Bamutaze, Y., Gidudu, A., Arinaitwe, G. and Voda, M., 2025. Dynamics and associations of selected agrometeorological variables in Robusta growing regions of Uganda. Agricultural Water Management, 307, p.109257.*|
|Health Care|*Ünal, E. and Yılmaz, S., 2025, March. Healthcare Sector Dynamics in Turkey (2002–2022): Trends, Breakpoints, and Policy Implications (Privatization in the Hospital Sector). In Healthcare (Vol. 13, No. 6, p. 622). MDPI.*|
|Rock Engineering|*Long, S., Yue, Z., Yue, W.V., Hu, H., Feng, Y., Yan, Y. and Xie, X., 2025. Identification of rock layer interface characteristics using drilling parameters. Rock Mechanics and Rock Engineering, 58(1), pp.1071-1098.*|
|Energy|*Leiria, D., Johra, H., Anoruo, J., Praulins, I., Piscitelli, M.S., Capozzoli, A., Marszal-Pomianowska, A. and Pomianowski, M.Z., 2025. Is it returning too hot? Time series segmentation and feature clustering of end-user substation faults in district heating systems. Applied Energy, 381, p.125122.*|
|Psychology|*Huang, X., Nguyen, A. and Lajoie, S.P., 2025. Examining socially shared regulation of learning in medical training: the interplay of heart rate change points on regulatory interactions: Huang et al. European Journal of Psychology of Education, 40(3), p.93.*|
|Paleoclimatology|*Anastasia Zhuravleva et al., 2023. Caribbean salinity anomalies contributed to variable North Atlantic circulation and climate during the Common Era.  Science Advances,   DOI:10.1126/sciadv.adg2639*|
|   Population Ecology  | *Henderson, P. A. (2021). [Southwood's Ecological Methods (5th edition)](https://www.google.com/books/edition/Southwood_s_Ecological_Methods/snEhEAAAQBAJ?hl=en&gbpv=1&pg=PA473&printsec=frontcover). Oxford University Press., page 475-476*|
|Climatology|*Webster, M.A., Riihelä, A., Kacimi, S., Ballinger, T.J., Blanchard-Wrigglesworth, E., Parker, C.L. and Boisvert, L., 2024. Summer snow on Arctic sea ice modulated by the Arctic Oscillation. Nature Geoscience, pp.1-8.*|
|Cardiology|*Ozier, D., Rafiq, T., de Souza, R. and Singh, S.M., 2023. [Use of Sacubitril/Valsartan Prior to Primary Prevention Implantable Cardioverter Defibrillator Implantation](https://doi.org/10.1016/j.cjco.2022.10.005). CJC Open.*|
|Spatial Ecology|*Laurin, G.V., Cotrina-Sanchez, A., Belelli-Marchesini, L., Tomelleri, E., Battipaglia, G., Cocozza, C., Niccoli, F., Kabala, J.P., Gianelle, D., Vescovo, L. and Da Ros, L., 2024. Comparing ground below-canopy and satellite spectral data for an improved and integrated forest phenology monitoring system. Ecological Indicators, 158, p.111328.*|
|Anthropocene Science|*Thomas, E.R., Vladimirova, D.O., Tetzner, D.R., Emanuelsson, D.B., Humby, J., Turner, S.D., Rose, N.L., Roberts, S.L., Gaca, P. and Cundy, A.B., 2023. The Palmer ice core as a candidate Global boundary Stratotype Section and Point for the Anthropocene series. The Anthropocene Review, p.20530196231155191.*|
|Biomedical Engineering|*Saghbiny, E., Da Silva, J., Leblanc, L., Bobbio, C., Morel, G.G., Vialle, R. and Tamadazte, B., 2023, September. Breach detection in spine surgery based on cutting torque with ex-vivo experimental validation. In Conference on New Technologies for Computer and Robot Assisted Surgery.*|
|Political Science|*Reuning, K., Whitesell, A. and Hannah, A.L., 2022. [Facebook algorithm changes may have amplified local republican parties](https://doi.org/10.1177/20531680221103809). Research & Politics, 9(2), p.20531680221103809.*|
|Food Science|*Zaytsev, V., Tutukina, M.N., Chetyrkina, M.R., Shelyakin, P.V., Ovchinnikov, G., Satybaldina, D., Kondrashov, V.A., Bandurist, M.S., Seilov, S., Gorin, D.A. and Fedorov, F.S., 2024. Monitoring of meat quality and change-point detection by a sensor array and profiling of bacterial communities. Analytica Chimica Acta, p.343022.*|
|Ecology|*Dashti, H., Chen, M., Smith, B., Zhao, K. and Moore, D., 2024. Rethinking Ecosystems Disturbance Recovery: what it was or what it could have been?. Geophysical Research Letters*|
|Spatial Hydrology|*Wang, Yiming, Xuesong Zhang, Kaiguang Zhao, and Debjani Singh. "Streamflow in the United States: Characteristics, trends, regime shifts, and extremes." Scientific Data 11, no. 1 (2024): 788.*|
|Business Science|*Li, Z. and Tian, Y., 2024. Skewed multifractal cross-correlation between price and volume during the COVID-19 pandemic: Evidence from China and European carbon markets. Applied Energy, 371, p.123716.*|
|Ecography|*Smith, M.M. and Pauli, J.N., 2024. Small but connected islands can maintain populations and genetic diversity under climate change. Ecography, p.e07119.*|
|Economics|*Sapkota, B.P., 2024. Analysis of Climate Policy and Monetary Policy Nexus in the Norwegian Context: A DSGE Approach (Master's thesis, Norwegian University of Life Sciences).*|
|Cognitive Science|*Prein, J.C., Maurits, L., Werwach, A., Haun, D.B. and Bohn, M., 2024. Variation in gaze following across the life span: A process‐level perspective. Developmental Science, p.e13546.*|
|Neuroscience|*Aqel, K., Wang, Z., Peng, Y.B. and Maia, P.D., 2024. Reconstructing rodent brain signals during euthanasia with eigensystem realization algorithm (ERA). Scientific Reports, 14(1), p.12261.*|
|Glaciology|*Ramón, C.L., Rueda, F.J., Priet‐Mahéo, M.C. and Andradóttir, H., 2024. The impact of deep glacial water diversions from a hydroelectric reservoir in the thermal dynamics of a sub-arctic lake. Journal of Hydrology, 635, p.131081.*|
|Quaternary Science|*Gibson, D.K., Bird, B.W., Finney, B.P. and Steinman, B.A., 2024. Holocene insolation and sea surface temperature influences on the polar front jet stream and precipitation in the midcontinental United States. Quaternary Science Reviews, 340, p.108865.*|
|Geography|*Lyu, R., Pang, J., Zhang, J. and Zhang, J., 2024. The impacts of disturbances on mountain ecosystem services: Insights from BEAST and Bayesian network. Applied Geography, 162, p.103143.*|
|Watershed Hydrology|*Sakizadeh, M., Milewski, A. and Sattari, M.T., 2023. Analysis of Long-Term Trend of Stream Flow and Interaction Effect of Land Use and Land Cover on Water Yield by SWAT Model and Statistical Learning in Part of Urmia Lake Basin, Northwest of Iran. Water, 15(4), p.690.*|
|Oceanography|*Oehlert, A.M., Hunter, H., Riopelle, C. and Purkis, S.J., 2023. Perturbation to North Atlantic Ocean‐Climate Dynamics Tripled Whitings Mud Production in the Bahamas. Journal of Geophysical Research: Oceans, 128(11), p.e2023JC020021.*|
| Hydraulic Engineering | *Xu, X., Yang, J., Ma, C., Qu, X., Chen, J. and Cheng, L., 2022. Segmented modeling method of dam displacement based on BEAST time series decomposition. Measurement, 202, p.111811.* |
|Social Media|*Barrie, C., Ketchley, N., Siegel, A. and Bagdouri, M., 2023. Measuring Media Freedom.*|
|Political Economy|*Benchimol, J. and Palumbo, L., 2023. Sanctions and Russian Online Prices.*|
|Physiology|*Shakeel, M., Brockmann, A. Temporal effects of sugar intake on fly local search and honey bee dance behaviour. J Comp Physiol A (2023). https://doi.org/10.1007/s00359-023-01670-6*|
|Injuries & Hazards|*Delavary, M., Kalantari, A.H., Mohammadzadeh Moghaddam, A., Fakoor, V., Lavallière, M. and Wilhelm Siebert, F., 2024. Road traffic mortality in Iran: longitudinal trend and seasonal analysis, March 2011-February 2020. International journal of injury control and safety promotion, 31(1), pp.125-137.*|
|Civil Engineering|*Langtry, M., Wichitwechkarn, V., Ward, R., Zhuang, C., Kreitmair, M.J., Makasis, N., Conti, Z.X. and Choudhary, R., 2024. Impact of data for forecasting on performance of model predictive control in buildings with smart energy storage. Energy and Buildings, p.114605.*|
|Ichthyology|*Kaeding, L.R., 2023. Climate-change and nonnative-piscivore impacts on a renowned Oncorhynchus metapopulation, requirements for metapopulation recovery, and similarities to anadromous salmonid metapopulations. Aquatic Sciences, 85(4), p.88.*|
|Remote Sensing|*Mulverhill, C., Coops, N.C. and Achim, A., 2023. Continuous monitoring and sub-annual change detection in high-latitude forests using Harmonized Landsat Sentinel-2 data. ISPRS Journal of Photogrammetry and Remote Sensing, 197, pp.309-319.*|
|Physical Chemistry|*Faran, M. and Bisker, G., 2023. Nonequilibrium Self-Assembly Time Forecasting by the Stochastic Landscape Method. The Journal of Physical Chemistry B.*|
|Biogeochemistry|*Dahl, M., Gullström, M., Bernabeu, I., Serrano, O., Leiva‐Dueñas, C., Linderholm, H.W., Asplund, M.E., Björk, M., Ou, T., Svensson, J.R. and Andrén, E., 2024. A 2,000‐year record of eelgrass (Zostera marina L.) colonization shows substantial gains in blue carbon storage and nutrient retention. Global Biogeochemical Cycles, 38(3), p.e2023GB008039.*|
|Mechanobiology|*Faran, M., Ray, D., Nag, S., Raucci, U., Parrinello, M. and Bisker, G., 2024. A Stochastic Landscape Approach for Protein Folding State Classification. Journal of Chemical Theory and Computation.*|
|Analytical Chemistry|*Simic, M., Neuper, C., Hohenester, U. and Hill, C., 2023. Optofluidic force induction as a process analytical technology. Analytical and Bioanalytical Chemistry, pp.1-11.*|
| Ecosystem Sciences | *Lyu, R., Zhao, W., Pang, J., Tian, X., Zhang, J. and Wang, N., 2022. [Towards a sustainable nature reserve management: Using Bayesian network to quantify the threat of disturbance to ecosystem services](https://doi.org/10.1016/j.ecoser.2022.101483). Ecosystem Services, 58, p.101483.* |
| Environmental Sciences| *Nickerson, S., Chen, G., Fearnside, P., Allan, C.J., Hu, T., de Carvalho, L.M. and Zhao, K., 2022. [Forest loss is significantly higher near clustered small dams than single large dams per megawatt of hydroelectricity installed in the Brazilian Amazon](https://iopscience.iop.org/article/10.1088/1748-9326/ac8236/meta). Environmental Research Letters.*|
|Geology| *Fan, X., Goeppert, N. and Goldscheider, N., 2023. Quantifying the historic and future response of karst spring discharge to climate variability and change at a snow-influenced temperate catchment in central Europe. Hydrogeology Journal, pp.1-17.*|
|Wildlife|*Smith, Matthew M., and Jonathan N. Pauli. "Connectivity maintains genetic diversity and population persistence within an archipelagic refugia even under declining lake ice." Mechanisms of species recovery for a forest carnivore in a changing landscape: 173.*|
|  Climate Sciences|*Duke, N.C., Mackenzie, J.R., Canning, A.D., Hutley, L.B., Bourke, A.J., Kovacs, J.M., Cormier, R., Staben, G., Lymburner, L. and Ai, E., 2022. [ENSO-driven extreme oscillations in mean sea level destabilise critical shoreline mangroves—An emerging threat](https://journals.plos.org/climate/article?id=10.1371/journal.pclm.0000037). PLOS Climate, 1(8), p.e000003*|
| Finance|*Candelaria, Christopher A., Shelby M. McNeill, and Kenneth A. Shores. (2022). What is a School Finance Reform? Uncovering the ubiquity and diversity of school finance reforms using a Bayesian changepoint estimator.(EdWorkingPaper: 22-587). Retrieved from Annenberg Institute at Brown University: https://doi.org/10.26300/4vey-3w10*|
| Public health|*Linnell, K., Fudolig, M., Schwartz, A., Ricketts, T.H., O'Neil-Dunne, J.P., Dodds, P.S. and Danforth, C.M., 2022. [Spatial changes in park visitation at the onset of the pandemic](https://journals.plos.org/globalpublichealth/article?id=10.1371/journal.pgph.0000766). arXiv preprint arXiv:2205.15937.*|
|Biometerology|*Li, Y., Liu, Y., Bohrer, G., Cai, Y., Wilson, A., Hu, T., Wang, Z. and Zhao, K., 2022. [Impacts of forest loss on local climate across the conterminous United States: Evidence from satellite time-series observations](https://doi.org/10.1016/j.scitotenv.2021.149651). Science of The Total Environment, 802, p.149651.*|
| Applied Math|*Ferguson, Daniel, and Francois G. Meyer. [Probability density estimation for sets of large graphs with respect to spectral information using stochastic block models](https://arxiv.org/abs/2207.02168). arXiv preprint arXiv:2207.02168 (2022).*|
|Transportation Science| *Delavary, M., Kalantari, A.H., Mohammadzadeh Moghaddam, A., Fakoor, V., Lavalliere, M. and Wilhelm Siebert, F., 2023. Road traffic mortality in Iran: longitudinal trend and seasonal analysis, March 2011-February 2020. International Journal of Injury Control and Safety Promotion, pp.1-12.*|
|Water quality|*He, Ziming, Jiayu Yao, Yancen Lu, and Danlu Guo. "Detecting and explaining long-term changes in river water quality in south eastern Australia." Hydrological Processes: e14741.*|
|Air quality|*Wu, S., Yao, J., Wang, Y. and Zhao, W., 2023. Influencing factors of PM2. 5 concentrations in the typical urban agglomerations in China based on wavelet perspective. Environmental Research, p.116641.*|
| Hydrology | *Zohaib, M. and Choi, M., 2020. [Satellite-based global-scale irrigation water use and its contemporary trends](https://doi.org/10.1016/j.rse.2022.113222). Science of The Total Environment, 714, p.136719.* |
| Energy Engineering |*Lindig, S., Theristis, M. and Moser, D., 2022. Best practices for photovoltaic performance loss rate calculations. Progress in Energy, 4(2), p.022003.*|
|Virology|*Shen, L., Sun, M., Song, S., Hu, Q., Wang, N., Ou, G., Guo, Z., Du, J., Shao, Z., Bai, Y. and Liu, K., 2022. The impact of anti-COVID19 nonpharmaceutical interventions on hand, foot, and mouth disease—A spatiotemporal perspective in Xi'an, northwestern China. Journal of medical virology.*|
| Pharmaceutical Sciences|*Patzkowski, M.S., Costantino, R.C., Kane, T.M., Nghiem, V.T., Kroma, R.B. and Highland, K.B., 2022. Military Health System Opioid, Tramadol, and Gabapentinoid Prescription Volumes Before and After a Defense Health Agency Policy Release. Clinical Drug Investigation, pp.1-8.*|
| Geography|*Cai, Y., Liu, S. and Lin, H., 2020. Monitoring the vegetation dynamics in the Dongting Lake Wetland from 2000 to 2019 using the BEAST algorithm based on dense Landsat time series. Applied Sciences, 10(12), p.4209.*|
| Oceanography|*Pitarch, J., Bellacicco, M., Marullo, S. and Van Der Woerd, H.J., 2021. [Global maps of Forel-Ule index, hue angle and Secchi disk depth derived from 21 years of monthly ESA Ocean Colour Climate Change Initiative data](https://essd.copernicus.org/articles/13/481/2021/). Earth System Science Data, 13(2), pp.481-490.*|
|Photovoltaics|*Micheli, L., Theristis, M., Livera, A., Stein, J.S., Georghiou, G.E., Muller, M., Almonacid, F. and Fernadez, E.F., 2021. Improved PV soiling extraction through the detection of cleanings and change points. IEEE Journal of Photovoltaics, 11(2), pp.519-526.*|
|Climate Sciences|*White, J.H., Walsh, J.E. and Thoman Jr, R.L., 2021. [Using Bayesian statistics to detect trends in Alaskan precipitation](https://doi.org/10.1002/joc.6946). International Journal of Climatology, 41(3), pp.2045-2059.*|
|Field Hydrology|*Merk, M., Goeppert, N. and Goldscheider, N., 2021. [Deep desiccation of soils observed by long-term high-resolution measurements on a large inclined lysimeter](https://hess.copernicus.org/articles/25/3519/2021/). Hydrology and Earth System Sciences, 25(6), pp.3519-3538.*|
|Sports Science|*Roccetti, M. and Piconi, F., Minutes played by Under 21 players in Serie A: a descriptive analytics study.*|
|Forest Ecology|*Moreno-Fernandez, D., Viana-Soto, A., Camarero, J.J., Zavala, M.A., Tijerin, J. and Garcia, M., 2021. Using spectral indices as early warning signals of forest dieback: The case of drought-prone Pinus pinaster forests. Science of The Total Environment, 793, p.148578.*|
|Petroleum  Engineering|*Pan, Y., Bi, R., Yang, S., Lyu, Z. and Ju, X., 2024, February. Application of a Bayesian Ensemble Algorithm for Automated Production Diagnostic of Gas Wells with Plunger-Lift. In International Petroleum Technology Conference (p. D011S029R008). IPTC.*|
|Atmospheric Sciences|*Tingwei, C., Tingxuan, H., Bing, M., Fei, G., Yanfang, X., Rongjie, L., Yi, M. and Jie, Z., 2021. Spatiotemporal pattern of aerosol types over the Bohai and Yellow Seas observed by CALIOP. Infrared and Laser Engineering, 50(6), p.20211030.*|
|Terrestrial ecology|*Dashti, H., Pandit, K., Glenn, N.F., Shinneman, D.J., Flerchinger, G.N., Hudak, A.T., de Graaf, M.A., Flores, A., Ustin, S., Ilangakoon, N. and Fellows, A.W., 2021. [Performance of the ecosystem demography model (EDv2. 2) in simulating gross primary production capacity and activity in a dryland study area](https://doi.org/10.1016/j.agrformet.2020.108270). Agricultural and Forest Meteorology, 297, p.108270.*|
|Statistics|*Storath, M. and Weinmann, A., 2023. Smoothing splines for discontinuous signals. Journal of Computational and Graphical Statistics, (just-accepted), pp.1-26.*|
|Environmental Engineering|*Bainbridge, R., Lim, M., Dunning, S., Winter, M.G., Diaz-Moreno, A., Martin, J., Torun, H., Sparkes, B., Khan, M.W. and Jin, N., 2022. Detection and forecasting of shallow landslides: lessons from a natural laboratory. Geomatics, Natural Hazards and Risk, 13(1), pp.686-704.*|
|Fishery|*Theis, S., Wallace, A.,, Poesch, M.,  Portiss, R. & Ruppert, J.. (2024). Balancing boat-electrofishing sampling effort against costs for nearshore fish communities in the Toronto waterfront, Lake Ontario. Fisheries Management and Ecology. 10.1111/fme.12733.*|
|Hydrology|*Yang, X., Tian, S., You, W. and Jiang, Z., 2021. Reconstruction of continuous GRACE/GRACE-FO terrestrial water storage anomalies based on time series decomposition. Journal of Hydrology, 603, p.127018.*|
|Landscape Ecology|*Adams, B.T., Matthews, S.N., Iverson, L.R., Prasad, A.M., Peters, M.P. and Zhao, K., 2021. Spring phenological variability promoted by topography and vegetation assembly processes in a temperate forest landscape. Agricultural and Forest Meteorology, 308, p.108578.*|

