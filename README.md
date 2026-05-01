# EVT-Portfolio-Risk
Title: EVT-Based Portfolio Risk & Stress Testing

 Overview
This project applies Extreme Value Theory (EVT) to model extreme financial losses and improve tail risk estimation.
It extends a previous stress testing framework by addressing limitations in traditional risk models.

 Objectives
•	Model extreme portfolio losses 
•	Estimate VaR and Expected Shortfall 
•	Compare GEV vs GPD approaches 
•	Enhance stress testing using EVT 

Key Features
•	Block maxima modeling (GEV) 
•	Peaks-over-threshold modeling (GPD) 
•	EVT-based VaR and Expected Shortfall 
•	Tail risk simulation 
•	Model validation and back testing 

Methodology
1.	Data preprocessing 
2.	Portfolio construction 
3.	GEV modeling 
4.	GPD modeling 
5.	Tail risk estimation 
6.	Stress testing extension 

 Key Results
•	Heavy-tailed behavior observed 
•	EVT improves extreme loss estimation 
•	GPD provides more stable tail modeling 
•	Traditional models underestimate risk 

 Project Evolution
•	Previous Project: Portfolio Stress Testing & Scenario Analysis 
o	Used historical simulation and Monte Carlo 
o	Identified limitations of normality assumptions 

•	Current Project: EVT Risk Modeling 
o	Models tail risk explicitly using GEV and GPD 
o	Improves estimation of extreme losses 

•	Next Project: FX Stress Testing (Planned) 
o	Apply EVT to USD/KES 
o	Incorporate macroeconomic shocks 
o	Model capital flow and currency risk

Repository Structure
EVT-Portfolio-Risk/
│── data/
│── scripts/
│── plots/
│── report/
│── README.md

 Tools
•	R 
•	dplyr, lubridate, EnvStats, QRM 

 Applications
•	Portfolio risk management 
•	Stress testing 
•	Tail risk modeling 
•	Financial engineering 

Key Insight
Extreme losses are not well captured by traditional models. EVT provides a more realistic framework for modeling financial risk.


Financial Engineering & Risk Analytics Project.

## Author
Maina Silvia
