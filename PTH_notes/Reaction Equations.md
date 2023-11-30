$\frac{dP}{dt}=k_{P,b}-P\cdot(s_{P,d}\cdot[Hb]+k_{R,d0})-2P\cdot ln(2)\cdot (t_{P,a,0}-\frac{t_{R,a,max}}{1+e^{-k_{R,a}\cdot(Hkt-Hkt_0)}})^{-1}$ 
$\frac{dR}{dt}=2P\cdot ln(2)\cdot (t_{P,a,0}-\frac{t_{R,a,max}}{1+e^{-k_{R,a}\cdot(Hkt-Hkt_0)}})^{-1}-2R\cdot ln(2)\cdot(\frac{t_{R,a,max}}{1+e^{-k_{R,a}\cdot(Hkt-Hkt_0)}})^{-1}-R\cdot M\cdot k_{inf}-R\cdot(k_{R,d}+s_{BH}\cdot J_{oiE}),k_{R,d}=0$
$\frac{dE}{dt}=R\cdot k_{R,a}-E\cdot M\cdot k_{inf}-E\cdot(k_{E,d}+s_{BH}\cdot oiE_{12} \cdot k_{oiE,d})$

**Infection**
$\frac{diE}{dt}= M\cdot(E+R)-iE\cdot (k_{iE,d,0} + (1-k_{iE,pit,frac}) \cdot \frac{k_{iE, ART_{max}} \cdot (ART^{h_{ART}}}{ART^{h_{ART}}+ID_{50}^{h_{ART}}})-iE\cdot k_{iE,rupture}$
$\frac{dM}{dt}=16\cdot iE\cdot k_{iE,rupture}-M\cdot k_{M,d}$
$\frac{doiE}{dt}=iE\cdot k_{iE,pit}-\sum\limits_{i} \frac{doiE_i}{dt}=>\frac{doiE_i}{dt}=oiE_i\cdot k_{oiE,d}$ or just $oiE_{12}\cdot k_{oiE,d}$ as influxes and outfluxes cancel another out?
**Medication**
$\frac{dART}{dt}=ART\cdot\frac{ln(2)}{t_{ART,decay}}$
#TODO: Flux for cell deaths which release LDH
$\frac{dLDH}{dt}=\frac{LDH_{Blood}}{V_{Blood}}\cdot \sum J_{d, of RBC cells}-LDH\cdot \frac{ln(2)}{t_{hl,LDH,d}}$, no idea of degradation time

## Parameters
**Erythopoesis**
- Precursors half life aging:  $t_{P,a}=t_{P,a,0}-\frac{t_{R,a,max}}{1+e^{-k_{R,a}\cdot(Hkt-Hkt_0)}}$ with $t_{P,a,0}=11$ #TODO why 11 not 12?
- Precursors death rate: #TODO Where are these values coming from?
	$s_{P,d} = 0.00071535$   slope of P death increase
	$k_{P,d0} =  0.48924947$ default death rate of Precursors
- Reticocyotes aging time: $t_{R,a}=t_{R,a,0}+\frac{t_{R,a,max}}{1+e^{-k_{R,a}\cdot(Hkt-Hkt_0)}}$, with $t_{R,a0}=0$
**Infection**
- $k_{iE,rupture}=\frac{2\cdot ln(2)}{t_{iE,rupture}}$ , $t_{iE,rupture=2\lor 4}$ CHECK/MAXIM: in 1/days, Austin(1998)
- $k_{iE,pit}=k_{iE,pit,0} + k_{iE,pit,frac} \cdot \frac{k_{iE, ART_{max}}\cdot ART^{h_{ART}}}{ART^{h_{ART}}+ID_{50}^{h_{ART}}}$
	- $k_{iE,pit,0}= 0$, assumption no pitting without ART
	- $k_{iE,pit,frac} = 0.33$ fraction of iE that are pitted, here assumption $\frac{1}{3}$, as only ring stage pitted survive to be oiE, trophozoite and schizont stage #TODO find citation for ring stage only pitting
- Bystander Haemolysis
	- $BH=s_{BH}\cdot J_{oiE,d}$  with $s_{BH}=k_{E,d}\cdot x$, with x as impact factor of dying oiE flux
	- Not implemented Idea: Hill Kinetics for BH as this can reach BH=0, as the linear function above, which logistic function cannot:
		- $BH=\frac{J_{oiE,d}^h}{BH_{max}\cdot J_{oiE,d}^h+K_m^h}$
- Linear chain trick used to model oiE (life span 7-20 days), as Florian talked about delayed death, with 12 species and $k_{oiE,d}\approx 0.963$ both values from fitting gamma distribution

**Medication**
- As a model event at t=0, ART is set to 400, since medication given to patient at the beginning of data collection
- $t_{ART,d}=1/12 days$ or 2h $\approx$ matches findings in N. J. White, [Clinical Pharmacokinetics and Pharmacodynamics of Artemisinin and Derivatives](https://doi.org/10.1016/0035-9203(94)90471-5), (June 1994) 
- $LDH_{Blood}=140-280 \frac{U}{L}$ only value of LDH,  [link](https://www.ncbi.nlm.nih.gov/books/NBK557536/?report=printable)
### Steady state
In model or estimation calculated during initialization of model .ant file for Steady State consideration of erythropoesis
$LDH_{Blood}=\frac{[LDH] \cdot \frac{ln(2)} {t_{hlf,LDH,d}} \cdot Vol_{blood}}{[E]\cdot \frac{2 \cdot ln(2)}{t_{E,d}} +\underbrace{J_{R,d}}_{=0}}$ with $t_{hl,LDH,d}$ = 4  bad citation for **in** erythrocytes [3-5 Tage](https://www.medicoconsult.de/ldh/) #TODO find better citation
Birth rate of precursor cells, calculated from steady state and then set for whole time 
$k_{P,b}= P\cdot(s_{P,d}\cdot Hb_{init}+k_{R,d0})-P(\frac{2\cdot ln(2)}{t_{P,a}})$ 
$Hb_{init} = \frac{(Vol_E \cdot E_{init} * Hb_{conc,E} + Vol_R \cdot R_{init} \cdot Hb_{conc,R})}{(10 \cdot Vol_{blood})}$ in model change of concentration unit from g/l to g/dl, thus div. by 10

Initial erythropoesis composition
$E_{init} = \frac{Hkt_{init} \cdot Vol_{blood}}{Vol_{E} + \frac{t_{R,a,init}}{t_{E,d}} \cdot Vol_{R}}$
$R_{init} = E_{init} \cdot \frac{t_{R,a,init}}{t_{E,d}}$
$P_{init} = \frac{1}{1024}\cdot R_{init}\cdot \frac{t_{P,a,init}} {t_{R,a,init}}$


