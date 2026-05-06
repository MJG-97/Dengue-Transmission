# SIR-SI Dengue Transmission Model with Climate Forcing
# Modeling the transmission of dengue fever in Mayote

## Description

This repository contains a compartmental SIR-SI mathematical model (5 compartments) for dengue fever transmission dynamics, incorporating climatic forcing to capture the impact of temperature and precipitation on the transmission rate. Three climate forcing models are compared: constant, Lambrechts (2011), and Gaussian.

## Code Structure

- Data loading and preparation
- Climate functions (Lambrechts, Gaussian with precipitation)
- SIR-SI ODE model (5 equations)
- Weekly incidence prediction functions
- Parameter estimation via Maximum Likelihood Estimation (Negative Binomial)
- Optimization and model comparison (AIC)
- Visualizations (observed/predicted incidences, climate-cases correlations, β(t) dynamics)

## Compared Models

| Model | Description |
|-------|-------------|
| Constant | β(t) constant over time |
| Lambrechts | β(t) depends on temperature via a function from the literature |
| Gaussian | β(t) depends on temperature (Gaussian kernel) and precipitation, with a lag of l weeks |

## Estimation Methods

- **MLE**: Maximum Likelihood Estimation with Negative Binomial distribution (default)
- **WLS**: Weighted Least Squares (available)

## Data

The data used for this project are confidential and not included in this repository. The code expects a CSV file with the following columns:

| Column | Description |
|--------|-------------|
| `n_cases` | Weekly case counts |
| `temp_mean_week` | Weekly mean temperature (°C) |
| `rain_mean_week` | Weekly mean precipitation (mm) |

## Prerequisites
Required R packages:
```r
dplyr, ggplot2, tidyr, deSolve, gridExtra, readr, lubridate, scales, knitr, numDeriv
```

## References
```bib

@article{liebig_forecasting_2021,
	title = {Forecasting the probability of local dengue outbreaks in {Queensland}, {Australia}},
	volume = {34},
	issn = {1755-4365},
	url = {https://www.sciencedirect.com/science/article/pii/S1755436520300426},
	doi = {10.1016/j.epidem.2020.100422},
	abstract = {The global incidence of dengue is increasing, and many previously unaffected areas have reported local cases of the vector-borne disease in recent years. For the effective containment of local outbreaks health authorities rely on the prompt notification of new cases. However, due to severe under-reporting and misdiagnosis, non-endemic countries face difficulties in containing local outbreaks, and the possibility of dengue becoming endemic. Outbreak control measures in non-endemic countries are largely reactive and health authorities would benefit from a universal early warning system that forecasts the probability of dengue outbreaks for given times and locations. We develop a model that establishes a link between pre- and post-border risk of dengue outbreaks. Specifically, we predict the probability of travellers importing dengue from other countries as well as the probability of those travellers causing local outbreaks. Our model can act as an early warning system, forecasting likely times and places of dengue outbreaks. We run our model for the Australian state of Queensland over a period of twelve years. Our results reveal the airports where dengue infected travellers are most likely to arrive and geographic locations associated with high outbreak probabilities. Our results can be used by health authorities to better utilise prevention and control resources and lead to the development of new prevention measures.},
	journal = {Epidemics},
	author = {Liebig, Jessica and de Hoog, Frank and Paini, Dean and Jurdak, Raja},
	month = mar,
	year = {2021},
	keywords = {Dengue, Early warning system, Forecasting disease spread, Human mobility},
	pages = {100422},
	file = {ScienceDirect Full Text PDF:C\:\\Users\\miche\\Zotero\\storage\\MKD68W6Y\\Liebig et al. - 2021 - Forecasting the probability of local dengue outbreaks in Queensland, Australia.pdf:application/pdf;ScienceDirect Snapshot:C\:\\Users\\miche\\Zotero\\storage\\6E2N77M4\\S1755436520300426.html:text/html},
}

@article{orellano_potential_2017,
	title = {Potential occurrence of {Zika} from subtropical to temperate {Argentina} considering the basic reproduction number ({R0})},
    journal = {Revista Panamericana de Salud Pública},
	author = {Orellano, Pablo and Vezzani, Darío and Quaranta, Nancy and Cionco, Rodolfo and Reynoso, Julieta and Salomon, Oscar},
    month = sep,
	year = {2017},
	pmid = {29466517},
	pmcid = {PMC6650626},
	pages = {e120},
	volume = {41},
	issn = {1020-4989},
	url = {https://pmc.ncbi.nlm.nih.gov/articles/PMC6650626/},
	doi = {10.26633/RPSP.2017.120},
	abstract = {Objective.
To assess the potential occurrence of Zika transmission throughout Argentina by the mosquito Aedes aegypti considering the basic reproduction number (R0).

Methods.
A model originally developed for dengue was adapted for Zika. R0 was estimatedas a function of seven parameters, three of them were considered temperature-dependent. Seasonal Zika occurrence was evaluated in 9 locations representing different climatic suitability for the vector. Data of diary temperatures were extracted and included in the model. A threshold of R0 = 1 was fixed for Zika occurrence. Sensitivity analyses were performed to evaluate the uncertainty around the results.

Results.
Zika transmission has the potential to occur in all studied locations at least in some moment of the year. In the northern region, transmission might be possible throughout the whole year or with an interruption in winter. The maximum R0 was estimated in 6.9, which means an average of 7 secondary cases from a primary case. The probabilistic sensitivity analysis showed that during winter the transmission can only be excluded in the southern fringe of geographic distribution of the vector and in part of central Argentina.

Conclusion.
Zika virus has the potential to be transmitted in Argentina throughout the current geographic range of the mosquito vector. Although the transmission would be mainly seasonal, the possibility of winter transmission cannot be excluded in northern and central Argentina, meaning that there is a potential endemic maintenance of the disease.},
	file = {Full Text PDF:C\:\\Users\\miche\\Zotero\\storage\\F6V6S8VN\\Orellano et al. - 2017 - Potential occurrence of Zika from subtropical to temperate Argentina considering the basic reproduct.pdf:application/pdf},
}


@article{mordecai_thermal_2019,
	title = {Thermal biology of mosquito-borne disease},
	volume = {22},
	copyright = {© 2019 The Authors Ecology Letters published by CNRS and John Wiley \& Sons Ltd},
	issn = {1461-0248},
	url = {https://onlinelibrary.wiley.com/doi/abs/10.1111/ele.13335},
	doi = {10.1111/ele.13335},
	abstract = {Mosquito-borne diseases cause a major burden of disease worldwide. The vital rates of these ectothermic vectors and parasites respond strongly and nonlinearly to temperature and therefore to climate change. Here, we review how trait-based approaches can synthesise and mechanistically predict the temperature dependence of transmission across vectors, pathogens, and environments. We present 11 pathogens transmitted by 15 different mosquito species – including globally important diseases like malaria, dengue, and Zika – synthesised from previously published studies. Transmission varied strongly and unimodally with temperature, peaking at 23–29ºC and declining to zero below 9–23ºC and above 32–38ºC. Different traits restricted transmission at low versus high temperatures, and temperature effects on transmission varied by both mosquito and parasite species. Temperate pathogens exhibit broader thermal ranges and cooler thermal minima and optima than tropical pathogens. Among tropical pathogens, malaria and Ross River virus had lower thermal optima (25–26ºC) while dengue and Zika viruses had the highest (29ºC) thermal optima. We expect warming to increase transmission below thermal optima but decrease transmission above optima. Key directions for future work include linking mechanistic models to field transmission, combining temperature effects with control measures, incorporating trait variation and temperature variation, and investigating climate adaptation and migration.},
	language = {en},
	number = {10},
	journal = {Ecology Letters},
	author = {Mordecai, Erin A. and Caldwell, Jamie M. and Grossman, Marissa K. and Lippi, Catherine A. and Johnson, Leah R. and Neira, Marco and Rohr, Jason R. and Ryan, Sadie J. and Savage, Van and Shocket, Marta S. and Sippy, Rachel and Stewart Ibarra, Anna M. and Thomas, Matthew B. and Villena, Oswaldo},
	year = {2019},
	keywords = {Arbovirus, climate change, malaria, mosquito, Ross River virus, temperature, thermal performance curve, West Nile virus, dengue virus, Zika virus},
	pages = {1690--1708},
	file = {Full Text PDF:C\:\\Users\\miche\\Zotero\\storage\\UNUZQKAE\\Mordecai et al. - 2019 - Thermal biology of mosquito-borne disease.pdf:application/pdf;Snapshot:C\:\\Users\\miche\\Zotero\\storage\\MJMYNP9M\\ele.html:text/html},
}

@article{andraud_dynamic_2012,
	title = {Dynamic {Epidemiological} {Models} for {Dengue} {Transmission}: {A} {Systematic} {Review} of {Structural} {Approaches}},
	volume = {7},
	issn = {1932-6203},
	shorttitle = {Dynamic {Epidemiological} {Models} for {Dengue} {Transmission}},
	url = {https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0049085},
	doi = {10.1371/journal.pone.0049085},
	abstract = {Dengue is a vector-borne disease recognized as the major arbovirose with four immunologically distant dengue serotypes coexisting in many endemic areas. Several mathematical models have been developed to understand the transmission dynamics of dengue, including the role of cross-reactive antibodies for the four different dengue serotypes. We aimed to review deterministic models of dengue transmission, in order to summarize the evolution of insights for, and provided by, such models, and to identify important characteristics for future model development. We identified relevant publications using PubMed and ISI Web of Knowledge, focusing on mathematical deterministic models of dengue transmission. Model assumptions were systematically extracted from each reviewed model structure, and were linked with their underlying epidemiological concepts. After defining common terms in vector-borne disease modelling, we generally categorised fourty-two published models of interest into single serotype and multiserotype models. The multi-serotype models assumed either vector-host or direct host-to-host transmission (ignoring the vector component). For each approach, we discussed the underlying structural and parameter assumptions, threshold behaviour and the projected impact of interventions. In view of the expected availability of dengue vaccines, modelling approaches will increasingly focus on the effectiveness and cost-effectiveness of vaccination options. For this purpose, the level of representation of the vector and host populations seems pivotal. Since vector-host transmission models would be required for projections of combined vaccination and vector control interventions, we advocate their use as most relevant to advice health policy in the future. The limited understanding of the factors which influence dengue transmission as well as limited data availability remain important concerns when applying dengue models to real-world decision problems.},
	language = {en},
	number = {11},
	journal = {PLOS ONE},
	author = {Andraud, Mathieu and Hens, Niel and Marais, Christiaan and Beutels, Philippe},
	month = nov,
	year = {2012},
	note = {Publisher: Public Library of Science},
	keywords = {Death rates, Dengue fever, Epidemiology, Immunity, Infectious disease control, Mosquitoes, Vaccination and immunization, Vector-borne diseases},
	pages = {e49085},
	file = {Full Text PDF:C\:\\Users\\miche\\Zotero\\storage\\FMDZV8WL\\Andraud et al. - 2012 - Dynamic Epidemiological Models for Dengue Transmission A Systematic Review of Structural Approaches.pdf:application/pdf},
}

@techreport{spfrance_dengue_mayotte_2020,
  author       = {{Santé publique France}},
  title        = {Le point épidémio : Dengue à Mayotte -- Point de situation au 10 décembre 2020},
  institution  = {Santé publique France},
  year         = {2020},
  date         = {2020-12-10},
  type         = {Bulletin épidémiologique},
  address      = {Mamoudzou, Mayotte, France},
  url          = {https://www.santepubliquefrance.fr/regions/ocean-indien/publications/#tabs}
}

@misc{ARS_ORS_2025_Mayotte,
  author       = {{Agence Régionale de la Santé (ARS)} and {Observatoire Régional de la Santé (ORS)}},
  title        = {Maladies vectorielles et comportements associés à Mayotte},
  year         = {2025},
  url          = {https://www.mayotte.ars.sante.fr/unono-ulanga-maladies-vectorielles-et-comportements}
}

@article{bhatt_current_2021,
	title = {Current {Understanding} of the {Pathogenesis} of {Dengue} {Virus} {Infection}},
	volume = {78},
	issn = {1432-0991},
	doi = {10.1007/s00284-020-02284-w},
	abstract = {The pathogenesis of dengue virus infection is attributed to complex interplay between virus, host genes and host immune response. Host factors such as antibody-dependent enhancement (ADE), memory cross-reactive T cells, anti-DENV NS1 antibodies, autoimmunity as well as genetic factors are major determinants of disease susceptibility. NS1 protein and anti-DENV NS1 antibodies were believed to be responsible for pathogenesis of severe dengue. The cytokine response of cross-reactive CD4+ T cells might be altered by the sequential infection with different DENV serotypes, leading to further elevation of pro-inflammatory cytokines contributing a detrimental immune response. Fcγ receptor-mediated antibody-dependent enhancement (ADE) results in release of cytokines from immune cells leading to vascular endothelial cell dysfunction and increased vascular permeability. Genomic variation of dengue virus and subgenomic flavivirus RNA (sfRNA) suppressing host immune response are viral determinants of disease severity. Dengue infection can lead to the generation of autoantibodies against DENV NS1antigen, DENV prM, and E proteins, which can cross-react with several self-antigens such as plasminogen, integrin, and platelet cells. Apart from viral factors, several host genetic factors and gene polymorphisms also have a role to play in pathogenesis of DENV infection. This review article highlights the various factors responsible for the pathogenesis of dengue and also highlights the recent advances in the field related to biomarkers which can be used in future for predicting severe disease outcome.},
	language = {eng},
	number = {1},
	journal = {Current Microbiology},
	author = {Bhatt, Puneet and Sabeena, Sasidharan Pillai and Varma, Muralidhar and Arunkumar, Govindakarnavar},
	month = jan,
	year = {2021},
	pmid = {33231723},
	pmcid = {PMC7815537},
	keywords = {Antibodies, Viral, Antibody-Dependent Enhancement, Dengue, Dengue Virus, Humans, Virus Diseases},
	pages = {17--32},
	file = {Texte intégral:C\:\\Users\\miche\\Zotero\\storage\\7QI57MJK\\Bhatt et al. - 2021 - Current Understanding of the Pathogenesis of Dengue Virus Infection.pdf:application/pdf},
}


@article{Andronico2024,
  title = {Comparing the Performance of Three Models Incorporating Weather Data to Forecast Dengue Epidemics in Reunion Island, 2018–2019},
  author = {Andronico, Alessio and Menudier, Luce and Salje, Henrik and Vincent, Muriel and Paireau, Juliette and de Valk, Henriette and Gallian, Pierre and Pastorino, Boris and Brady, Oliver and de Lamballerie, Xavier and Lazarus, Clément and Paty, Marie-Claire and Vilain, Pascal and Noel, Harold and Cauchemez, Simon},
  journal = {The Journal of Infectious Diseases},
  volume = {229},
  number = {1},
  pages = {10--18},
  year = {2024},
  month = {January},
  publisher = {Oxford University Press},
  doi = {10.1093/infdis/jiad468},
  url = {https://academic.oup.com/jid/article/229/1/10/7438952}
}


@article{Chan2010,
  author    = {Chan, Miranda and Johansson, Michael A.},
  title     = {The Incubation Periods of Dengue Viruses},
  journal   = {PLOS ONE},
  year      = {2012},
  volume    = {7},
  number    = {11},
  pages     = {e50972},
  doi       = {10.1371/journal.pone.0050972},
  url       = {https://doi.org/10.1371/journal.pone.0050972},
  publisher = {Public Library of Science}
}

@online{INRS2025,
  author  = {Institut National de Recherche et de Sécurité (INRS)},
  title   = {Données épidémiologiques sur la dengue},
  year    = {2025},
  url     = {https://www.inrs.fr/publications/bdd/eficatt/fiche.html?refINRS=EFICATT_Dengue}
}

@article{Lambrechts2011,
  author    = {Lambrechts, L. and Paaijmans, K. P. and Fansiri, T. and Carrington, L. B. and Kramer, L. D. and Thomas, M. B. and Scott, T. W.},
  title     = {Impact of daily temperature fluctuations on dengue virus transmission by Aedes aegypti},
  journal   = {Proceedings of the National Academy of Sciences of the United States of America},
  volume    = {108},
  number    = {18},
  pages     = {7460--7465},
  year      = {2011},
  doi       = {10.1073/pnas.1101377108}
}

@article{Lee2018,
  title={Potential effects of climate change on dengue transmission dynamics in Korea},
  author={Lee, H and Kim, JE and Lee, S and Lee, CH},
  journal={PLOS ONE},
  volume={13},
  number={6},
  year={2018},
  doi={https://doi.org/10.1371/journal.pone.0199205}
}


@unpublished{manou-abi_spatio-temporal_2023,
	title = {Spatio-temporal modeling and risk ratio assessment of adult and egg mosquitoes abundance with consideration of environmental data in the island of {Mayotte}},
	url = {https://hal.science/hal-04224216},
	abstract = {This work provides a spatio-temporal statistical modeling for egg and adult Aedes mosquitoes count data with consideration of environmental data in the context of mosquito epidemiology. For a given spatio-temporal incomplete mosquito count data, we derive predictions assuming that all of the spatio-temporal dependence can be accounted by potential factors influencing the development of Aedes mosquitoes such as, rainfall, temperature (including many delay) and waste data. In this paper, after a data analysis with the entomological and environmental data, we apply a LASSO regression to perform a variable selection strategy (we are in the presence of a large number of explanatory variables). We highlight the relevant factors that explain the abundance of our egg and adult mosquitoes count data. We define a spatio-temporal risk ratio which is a probability of exceeding a given threshold value of mosquito abundance. We propose two spatio-temporal modeling approaches for the Aedes mosquitoes' count data. The first is based on a spatio-temporal kernel smoother and the second on a generalized additive model. The paper conclude with a detailed discussion that follows not only, the spatio-temporal prediction and model performance measures, but also the obtained spatio-temporal risk ratio in order to highlight potential space-time areas of threshold exceedances of mosquitoes' abundance where health services could apply vector surveillance and control measures. We also discuss a forthcoming work concerning the theoretical developments of a spatio-temporal model for simulation purposes and an R Shiny application called May'Aedes, that aims to be a flexible and efficient tool to predict spatio-temporal risk ratio.},
	author = {Manou-Abi, Solym and Logamou Seknewna, Lema and Dabo-Niang, Sophie and Balicchi, Julien and Idaroussi, Ambdoul-Bar and Saindou, Maoulide},
	month = oct,
	year = {2023},
	annote = {working paper or preprint},
	file = {HAL PDF Full Text:C\:\\Users\\miche\\Zotero\\storage\\IMNPFFZW\\Manou-Abi et al. - 2023 - Spatio-temporal modeling and risk ratio assessment of adult and egg mosquitoes abundance with consid.pdf:application/pdf},
}



@article{Morin2013,
  author    = {Cory W. Morin and Andrew C. Comrie},
  title     = {Regional and seasonal response of a West Nile virus vector to climate change},
  journal   = {Proceedings of the National Academy of Sciences},
  volume    = {110},
  number    = {39},
  pages     = {15620--15625},
  year      = {2013},
  publisher = {National Academy of Sciences},
  doi       = {10.1073/pnas.1307135110},
  url       = {https://www.pnas.org/doi/10.1073/pnas.1307135110}
}

@article{Nepomuceno2020,
  title={Incorporating vector control interventions in dengue transmission models: A methodological review},
  author={Nepomuceno, João and others},
  journal={Infectious Disease Modelling},
  volume={5},
  pages={75--89},
  year={2020},
  doi={10.1016/j.idm.2019.12.002}
}

@article{Newton1992,
  author = {Newton, Elizabeth A. and Reiter, Paul},
  title = {A model of the transmission of dengue fever with an evaluation of the impact of ultra-low volume (ULV) insecticide applications on dengue epidemics},
  journal = {American Journal of Tropical Medicine and Hygiene},
  volume = {47},
  number = {6},
  pages = {709--720},
  year = {1992},
  doi = {10.4269/ajtmh.1992.47.709}
}

@article{Perkins2015,
  author = {Perkins, T. A. and Metcalf, C. J. E. and Grenfell, B. T. and Tatem, A. J.},
  title = {Estimating Drivers of Autochthonous Transmission of Chikungunya Virus in its Invasion of the Americas},
  journal = {PLOS Currents Outbreaks},
  volume = {Edition 1},
  year = {2015},
  month = {February},
  day = {10},
  doi = {10.1371/currents.outbreaks.a4c7b6ac10e0420b1788c9767946d1fc},
  url = {https://doi.org/10.1371/currents.outbreaks.a4c7b6ac10e0420b1788c9767946d1fc}
}



@misc{OWID_Dengue_Data_2025,
  author       = {Esteban Ortiz-Ospina and Max Roser},
  title        = {Data Page: Dengue fever infections},
  howpublished = {Part of the publication \emph{Global Health}},
  year         = {2016},
  note         = {Data adapted from IHME, Global Burden of Disease. Archived online on October 1, 2025.},
  url          = {https://archive.ourworldindata.org/20251001-144546/grapher/dengue-incidence.html},
}

@article{RamirezSoto2023,
  title={SIR--SI model with a Gaussian transmission rate: Understanding the dynamics of dengue outbreaks in Lima, Peru},
  author={Ramírez-Soto, Max Carlos and others},
  journal={PLOS ONE},
  volume={18},
  number={4},
  pages={1--22},
  year={2023},
  doi={10.1371/journal.pone.0284263}
}

@article{Rui2024,
  title={A six-step framework for developing infectious disease models},
  author={Rui, Song and others},
  journal={Epidemics},
  volume={46},
  pages={100676},
  year={2024},
  doi={10.1016/j.epidem.2024.100676}
}

@article{SanchezVargas2018,
  author    = {Irma Sánchez-Vargas and Laura C. Harrington and Jeffrey B. Doty and William C. Black IV and Ken E. Olson},
  title     = {Demonstration of efficient vertical and venereal transmission of dengue virus type-2 in a genetically diverse laboratory strain of Aedes aegypti},
  journal   = {PLOS Neglected Tropical Diseases},
  volume    = {12},
  number    = {8},
  pages     = {e0006754},
  year      = {2018},
  doi       = {10.1371/journal.pntd.0006754},
  url       = {https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0006754}
}

@incollection{Tomashek2011,
  author    = {Kim M. Tomashek and Harold S. Margolis},
  title     = {Dengue: Epidemiology, Clinical Features, and Diagnosis},
  booktitle = {Dengue},
  editor    = {Scott B. Halstead},
  series    = {Tropical Medicine: Science and Practice},
  volume    = {5},
  pages     = {25--44},
  year      = {2011},
  publisher = {Imperial College Press},
  doi       = {10.1142/9781848165287_0002}
}

@misc{WHO2020,
  author    = {{World Health Organization}},
  title     = {Dengue and severe dengue},
  year      = {2020},
  url       = {https://www.who.int/news-room/fact-sheets/detail/dengue-and-severe-dengue}}


@article{esteva_analysis_1998,
	title = {Analysis of a dengue disease transmission model},
	volume = {150},
	issn = {0025-5564},
	url = {https://www.sciencedirect.com/science/article/pii/S0025556498100032},
	doi = {10.1016/S0025-5564(98)10003-2},
	abstract = {A model for the transmission of dengue fever in a constant human population and variable vector population is discussed. A complete global analysis is given, which uses the results of the theory of competitive systems and stability of periodic orbits, to establish the global stability of the endemic equilibrium. The control measures of the vector population are discussed in terms of the threshold condition, which governs the existence and stability of the endemic equilibrium.},
	number = {2},
	journal = {Mathematical Biosciences},
	author = {Esteva, Lourdes and Vargas, Cristobal},
	month = jun,
	year = {1998},
	keywords = {Competitive systems, Dengue disease, Endemic equilibrium, Global stability, Threshold number},
	pages = {131--151},
	file = {ScienceDirect Snapshot:C\:\\Users\\miche\\Zotero\\storage\\2VQGCYVJ\\S0025556498100032.html:text/html},
}


@article{focks_simulation_1995,
	title = {A {Simulation} {Model} of the {Epidemiology} of {Urban} {Dengue} {Fever}: {Literature} {Analysis}, {Model} {Development}, {Preliminary} {Validation}, and {Samples} of {Simulation} {Results}},
	volume = {53},
	shorttitle = {A {Simulation} {Model} of the {Epidemiology} of {Urban} {Dengue} {Fever}},
	url = {https://www.ajtmh.org/view/journals/tpmd/53/5/article-p489.xml},
	doi = {10.4269/ajtmh.1995.53.489},
	abstract = {"A Simulation Model of the Epidemiology of Urban Dengue Fever: Literature Analysis, Model Development, Preliminary Validation, and Samples of Simulation Results" published on Nov 1995 by The American Society of Tropical Medicine and Hygiene.},
	language = {EN},
	number = {5},
	journal = {The American Journal of Tropical Medicine and Hygiene},
	author = {Focks, Dana A. and Daniels, Eric and Haile, Dan G. and Keesling, James E.},
	month = nov,
	year = {1995},
	note = {Publisher: The American Society of Tropical Medicine and Hygiene
Section: The American Journal of Tropical Medicine and Hygiene},
	pages = {489--506},
}


@article{favier_influence_2005,
	title = {Influence of spatial heterogeneity on an emerging infectious disease: the case of dengue epidemics},
	volume = {272},
	issn = {0962-8452},
	shorttitle = {Influence of spatial heterogeneity on an emerging infectious disease},
	url = {https://doi.org/10.1098/rspb.2004.3020},
	doi = {10.1098/rspb.2004.3020},
	abstract = {The importance of spatial heterogeneity and spatial scales (at a village or neighbourhood scale) has been explored with individual-based models. Our reasoning is based on the Chilean Easter Island (EI) case, where a first dengue epidemic occurred in 2002 among the relatively small population localized in one village. Even in this simple situation, the real epidemic is not consistent with homogeneous models. Conversely, including contact heterogeneity on different scales (intra-households, inter-house, inter-areas) allows the recovery of not only the EI epidemiological curve but also the qualitative patterns of Brazilian urban dengue epidemic in more complex situations.},
	number = {1568},
	journal = {Proceedings of the Royal Society B: Biological Sciences},
	author = {Favier, Charly and Schmit, Delphine and Müller-Graf, Christine D.M and Cazelles, Bernard and Degallier, Nicolas and Mondet, Bernard and Dubois, Marc A},
	month = jun,
	year = {2005},
	pages = {1171--1177},
	file = {Snapshot:C\:\\Users\\miche\\Zotero\\storage\\XAQN9EJF\\rspb.2004.html:text/html},
}


@article{chowell_influence_2011,
	title = {The influence of geographic and climate factors on the timing of dengue epidemics in {Perú}, 1994-2008},
	volume = {11},
	issn = {1471-2334},
	url = {https://doi.org/10.1186/1471-2334-11-164},
	doi = {10.1186/1471-2334-11-164},
	abstract = {Dengue fever is a mosquito-borne disease that affects between 50 and 100 million people each year. Increasing our understanding of the heterogeneous transmission patterns of dengue at different spatial scales could have considerable public health value by guiding intervention strategies.},
	language = {en},
	number = {1},
	journal = {BMC Infectious Diseases},
	author = {Chowell, Gerardo and Cazelles, Bernard and Broutin, Hélène and Munayco, Cesar V.},
	month = jun,
	year = {2011},
	keywords = {climatic factors, community size, Dengue, dynamics, epidemic timing, Perú, wavelet analysis, wavelet coherence},
	pages = {164},
	file = {Full Text PDF:C\:\\Users\\miche\\Zotero\\storage\\DCETZWLK\\Chowell et al. - 2011 - The influence of geographic and climate factors on the timing of dengue epidemics in Perú, 1994-2008.pdf:application/pdf},
}


@article{picinini_freitas_zika_2024,
	title = {Zika emergence, persistence, and transmission rate in {Colombia}: a nationwide application of a space-time {Markov} switching model},
	volume = {14},
	issn = {2045-2322},
	shorttitle = {Zika emergence, persistence, and transmission rate in {Colombia}},
	doi = {10.1038/s41598-024-59976-7},
	abstract = {Zika, a viral disease transmitted to humans by Aedes mosquitoes, emerged in the Americas in 2015, causing large-scale epidemics. Colombia alone reported over 72,000 Zika cases between 2015 and 2016. Using national surveillance data from 1121 municipalities over 70 weeks, we identified sociodemographic and environmental factors associated with Zika's emergence, re-emergence, persistence, and transmission intensity in Colombia. We fitted a zero-state Markov-switching model under the Bayesian framework, assuming Zika switched between periods of presence and absence according to spatially and temporally varying probabilities of emergence/re-emergence (from absence to presence) and persistence (from presence to presence). These probabilities were assumed to follow a series of mixed multiple logistic regressions. When Zika was present, assuming that the cases follow a negative binomial distribution, we estimated the transmission intensity rate. Our results indicate that Zika emerged/re-emerged sooner and that transmission was intensified in municipalities that were more densely populated, at lower altitudes and/or with less vegetation cover. Warmer temperatures and less weekly-accumulated rain were also associated with Zika emergence. Zika cases persisted for longer in more densely populated areas with more cases reported in the previous week. Overall, population density, elevation, and temperature were identified as the main contributors to the first Zika epidemic in Colombia. We also estimated the probability of Zika presence by municipality and week, and the results suggest that the disease circulated undetected by the surveillance system on many occasions. Our results offer insights into priority areas for public health interventions against emerging and re-emerging Aedes-borne diseases.},
	language = {eng},
	number = {1},
	journal = {Scientific Reports},
	author = {Picinini Freitas, Laís and Douwes-Schultz, Dirk and Schmidt, Alexandra M. and Ávila Monsalve, Brayan and Salazar Flórez, Jorge Emilio and García-Balaguera, César and Restrepo, Berta N. and Jaramillo-Ramirez, Gloria I. and Carabali, Mabel and Zinszer, Kate},
	month = may,
	year = {2024},
	pmid = {38693192},
	pmcid = {PMC11063144},
	keywords = {Aedes, Animals, Bayes Theorem, Colombia, Disease Outbreaks, Humans, Markov Chains, Mosquito Vectors, Zika Virus, Zika Virus Infection},
	pages = {10003},
	file = {Texte intégral:C\:\\Users\\miche\\Zotero\\storage\\69MEQJS8\\Picinini Freitas et al. - 2024 - Zika emergence, persistence, and transmission rate in Colombia a nationwide application of a space-.pdf:application/pdf},
}

@article{chowell_effective_2016,
    title={The effective reproduction number R_t for infectious diseases},
    author={Chowell, Gerardo and others},
    journal={Handbook of Statistics},
    volume={36},
    pages={87-117},
    year={2016}
}


@article{diekmann_construction_2010,
	title = {The construction of next-generation matrices for compartmental epidemic models},
	volume = {7},
	issn = {1742-5689},
	url = {https://pmc.ncbi.nlm.nih.gov/articles/PMC2871801/},
	doi = {10.1098/rsif.2009.0386},
	abstract = {The basic reproduction number ℛ0 is arguably the most important quantity in infectious disease epidemiology. The next-generation matrix (NGM) is the natural basis for the definition and calculation of ℛ0 where finitely many different categories of individuals are recognized. We clear up confusion that has been around in the literature concerning the construction of this matrix, specifically for the most frequently used so-called compartmental models. We present a detailed easy recipe for the construction of the NGM from basic ingredients derived directly from the specifications of the model. We show that two related matrices exist which we define to be the NGM with large domain and the NGM with small domain. The three matrices together reflect the range of possibilities encountered in the literature for the characterization of ℛ0. We show how they are connected and how their construction follows from the basic model ingredients, and establish that they have the same non-zero eigenvalues, the largest of which is the basic reproduction number ℛ0. Although we present formal recipes based on linear algebra, we encourage the construction of the NGM by way of direct epidemiological reasoning, using the clear interpretation of the elements of the NGM and of the model ingredients. We present a selection of examples as a practical guide to our methods. In the appendix we present an elementary but complete proof that ℛ0 defined as the dominant eigenvalue of the NGM for compartmental systems and the Malthusian parameter r, the real-time exponential growth rate in the early phase of an outbreak, are connected by the properties that ℛ0 {\textgreater} 1 if and only if r {\textgreater} 0, and ℛ0 = 1 if and only if r = 0.},
	number = {47},
	urldate = {2026-04-28},
	journal = {Journal of the Royal Society Interface},
	author = {Diekmann, O. and Heesterbeek, J. A. P. and Roberts, M. G.},
	month = jun,
	year = {2010},
	pmid = {19892718},
	pmcid = {PMC2871801},
	pages = {873--885},
	file = {Full Text PDF:C\:\\Users\\miche\\Zotero\\storage\\BREZUVF8\\Diekmann et al. - 2010 - The construction of next-generation matrices for compartmental epidemic models.pdf:application/pdf},
}


@misc{bacaer_seuil_2006,
	title = {Le seuil épidémique des maladies à vecteurs avec saisonnalité},
	url = {https://hal.science/hal-01285109},
	abstract = {La leishmaniose cutanée est une maladie à vecteurs transmise aux humains par des phlébotomes. On développe dans cet article un modèle mathématique qui tient compte de la saisonnalité de la population de vecteurs et de la distribution de la période de latence entre l'infection et les symptômes chez les humains. On ajuste les paramètres à des données de la province de Chichaoua au Maroc. On propose aussi une généralisation de la définition de la reproductivité nette R0 adaptée aux environnements périodiques. On estime ce R0 numériquement pour l'épidémie à Chichaoua: R0 ≃ 1,94. Le modèle suggère que l'épidémie s'arrêterait si la population de vecteurs était divisée par (R0){\textasciicircum}2 ≃ 3,76.},
	urldate = {2026-04-28},
	author = {Bacaër, Nicolas and Guernaoui, Souad},
	year = {2006},
	doi = {10.1007/s00285-006-0015-0},
	note = {Pages: 421 - 436
Volume: 53},
	file = {HAL PDF Full Text:C\:\\Users\\miche\\Zotero\\storage\\R8RI7MF3\\Bacaër et Guernaoui - 2006 - Le seuil épidémique des maladies à vecteurs avec saisonnalité.pdf:application/pdf},
}

@article{andronico_comparing_2024,
	title = {Comparing the {Performance} of {Three} {Models} {Incorporating} {Weather} {Data} to {Forecast} {Dengue} {Epidemics} in {Reunion} {Island}, 2018–2019},
	volume = {229},
	issn = {0022-1899},
	url = {https://doi.org/10.1093/infdis/jiad468},
	doi = {10.1093/infdis/jiad468},
	abstract = {We developed mathematical models to analyze a large dengue virus (DENV) epidemic in Reunion Island in 2018–2019. Our models captured major drivers of uncertainty including the complex relationship between climate and DENV transmission, temperature trends, and underreporting. Early assessment correctly concluded that persistence of DENV transmission during the austral winter 2018 was likely and that the second epidemic wave would be larger than the first one. From November 2018, the detection probability was estimated at 10\%–20\% and, for this range of values, our projections were found to be remarkably accurate. Overall, we estimated that 8\% and 18\% of the population were infected during the first and second wave, respectively. Out of the 3 models considered, the best-fitting one was calibrated to laboratory entomological data, and accounted for temperature but not precipitation. This study showcases the contribution of modeling to strengthen risk assessments and planning of national and local authorities.},
	number = {1},
	urldate = {2026-04-30},
	journal = {The Journal of Infectious Diseases},
	author = {Andronico, Alessio and Menudier, Luce and Salje, Henrik and Vincent, Muriel and Paireau, Juliette and de Valk, Henriette and Gallian, Pierre and Pastorino, Boris and Brady, Oliver and de Lamballerie, Xavier and Lazarus, Clément and Paty, Marie-Claire and Vilain, Pascal and Noel, Harold and Cauchemez, Simon},
	month = jan,
	year = {2024},
	pages = {10--18},
	file = {Full Text PDF:C\:\\Users\\miche\\Zotero\\storage\\3GFS9FN4\\Andronico et al. - 2024 - Comparing the Performance of Three Models Incorporating Weather Data to Forecast Dengue Epidemics in.pdf:application/pdf;Snapshot:C\:\\Users\\miche\\Zotero\\storage\\WY5K4K5F\\jiad468.html:text/html},
}

