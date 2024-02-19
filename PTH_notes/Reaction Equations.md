$\frac{dP}{dt}=k_{P,b}-P\cdot \frac{( 1 + k_{P,d} * Hb^{r_{P,d}})} {a_{P,d}}-2P\cdot ln(2)\cdot (t_{P,a,0}-\frac{t_{R,a,max}}{1+e^{-k_{R,a}\cdot(Hkt-Hkt_0)}})^{-1}$ 
$\frac{dR}{dt}=2P\cdot ln(2)\cdot (t_{P,a,0}-\frac{t_{R,a,max}}{1+e^{-k_{R,a}\cdot(Hkt-Hkt_0)}})^{-1}-2R\cdot ln(2)\cdot(\frac{t_{R,a,max}}{1+e^{-k_{R,a}\cdot(Hkt-Hkt_0)}})^{-1}-R\cdot M\cdot k_{inf,R}-R\cdot(k_{R,d}+s_{BH}\cdot oiE_{12}),with k_{R,d}=0$
$\frac{dE}{dt}=2R\cdot ln(2)\cdot(\frac{t_{R,a,max}}{1+e^{-k_{R,a}\cdot(Hkt-Hkt_0)}})^{-1}-E\cdot M\cdot k_{inf}-E\cdot(k_{E,d}+s_{BH}\cdot oiE_{12})$
$\frac{dLDH}{dt}=\frac{LDH_{Blood}}{V_{Blood}}\cdot((P \cdot \frac{( 1 + k_{P,d} * Hb^{r_{P,d}})} {a_{P,d}})+R\cdot(k_{R,d}+s_{BH}\cdot oiE_{12})+E\cdot(k_{E,d}+s_{BH}\cdot oiE_{12} ))-LDH\cdot \frac{ln(2)}{t_{hl,LDH,d}}$, no real idea of degradation time of LDH

**Infection**
$\frac{diE}{dt}= M\cdot(E+R)-iE\cdot (k_{iE,d,0} + (1-k_{iE,pit,frac}) \cdot \frac{k_{iE, ART_{max}} \cdot (ART^{h_{ART}}}{ART^{h_{ART}}+ID_{50}^{h_{ART}}})-iE\cdot k_{iE,rupture}$
$\frac{dM}{dt}=16\cdot iE\cdot k_{iE,rupture}-M\cdot k_{M,d}$
$\frac{doiE}{dt}=iE\cdot k_{iE,pit}-oiE_{12}\cdot k_{oiE,d}$ , as only real death rest of oiE species age
**Medication**
$\frac{dART}{dt}=ART\cdot\frac{ln(2)}{t_{ART,decay}}$

## Parameters
**Erythropoiesis**
- Precursors half life aging:
	- $t_{P,a}=t_{P,a,0}-\frac{t_{R,a,max}}{1+e^{-k_{R,a}\cdot(Hkt-Hkt_0)}}$ with $t_{P,a,0}=7$ #TODO why 11 not 12
- Precursors death rate: 
	- hill function first mentioned in [here](https://pubmed.ncbi.nlm.nih.gov/7606142/): explanation usual enzymatic function and data linear in semi-logarithmic scale so exponential suspected
- Reticulocytes aging time:
	- $t_{R,a}=t_{R,a,0}+\frac{t_{R,a,max}}{1+e^{-k_{R,a}\cdot(Hkt-Hkt_0)}}$, with $t_{R,a0}=0$
-  $LDH_{Blood}=140-280 \frac{U}{L}$ only value of LDH,  [link](https://www.ncbi.nlm.nih.gov/books/NBK557536/?report=printable)
- $t_{hl,LDH,d}$ = 4  bad citation for **in** erythrocytes [3-5 Tage](https://www.medicoconsult.de/ldh/) #TODO find better citation

**Infection:**
- $k_{inf}:=$ infection rate, to be estimated and since #tropism $k_{inf,R}= trop \cdot k_{inf}$[# Erythrocyte tropism of malarial parasites: The reticulocyte appeal, Leong et al. 2022, review]([https://doi.org/10.3389/fmicb.2022.1022828](https://doi.org/10.3389/fmicb.2022.1022828))
- $k_{iE,rupture}=\frac{2\cdot ln(2)}{t_{iE,rupture}}$ , 
	- $t_{iE,rupture}=2$: in 1/days, anything higher than 2 makes little sense as it is used like the half life time meaning 1/2 maximal value reached at this time point but most observed rupture is synchronize and all rupture after 2 days. 
- $k_{iE,pit}=k_{iE,pit,0} + k_{iE,pit,frac} \cdot \frac{k_{iE, ART_{max}}\cdot ART^{h_{ART}}}{ART^{h_{ART}}+ID_{50}^{h_{ART}}}$
	- $k_{iE,pit,0}= 0$, assumption no pitting without ART, actually a little pitting happens
		- Pierre A. Buffet et al., [Retention of Erythrocytes in the Spleen: A Double-Edged Process in Human Malaria, Current Opinion in Hematology](https://doi.org/10.1097/MOH.0b013e32832a1d4b) (May 2009)
	- $k_{iE,pit,frac} = 0.33$ fraction of iE that are pitted, here assumption $\frac{1}{3}$, as only ring stage pitted iE survive to be oiE, 
		- #Problem random guess of number as only rings survive pitting but also only rings go there, so no relation to stage distribution of iE
		- The other stages for P.fal. are not as deformable but for P.vivax they are, so species specific pitting 
		- Pierre A. Buffet et al., [Retention of Erythrocytes in the Spleen: A Double-Edged Process in Human Malaria, Current Opinion in Hematology](https://doi.org/10.1097/MOH.0b013e32832a1d4b) (May 2009)
- $k_{M,d}=\frac{ln(2)}{t_{hlf,M}}$, 
	- with $t_{hlf,M}=30min$ estimated(guessed?) in Austin et al., [The Dynamics of Drug Action on the Within-Host Population Growth of Infectious Agents: Melding Pharmacokinetics with Pathogen Population Dynamics]( https://doi.org/10.1006/jtbi.1997.0438)(1998) 
	- other publication uses P.knowlesi only 10min $=t_{hlf}$ Johnson, J G et al. [“Factors affecting the ability of isolated Plasmodium knowlesi merozoites to attach to and invade erythrocytes.”]( doi:10.1017/s0031182000000998) _Parasitology_ (1980)
- Bystander Haemolysis
	- $BH=s_{BH}\cdot J_{oiE,d}$  with $s_{BH}=x\cdot k_{E,d}$, with x as impact factor of dying oiE flux
	- Not implemented Idea: Hill Kinetics for BH as this can reach BH=0, as the linear function above, which logistic function cannot:
		- $BH=\frac{J_{oiE,d}^h}{BH_{max}\cdot J_{oiE,d}^h+K_m^h}$
- Linear chain trick used to model oiE (life span 7-20 days), as Florian talked about delayed death, with 12 species and $k_{oiE,d}\approx 0.963$ both values from fitting gamma distribution

**Medication**:
- As a model event at t=0, ART is set to 400, since medication given to patient at the beginning of data collection
- $t_{ART,d}=1/12 days$ or 2h $\approx$ matches findings in N. J. White, [Clinical Pharmacokinetics and Pharmacodynamics of Artemisinin and Derivatives](https://doi.org/10.1016/0035-9203(94)90471-5), (June 1994) 
### Steady state
In model or estimation calculated during initialization of model .ant file for Steady State consideration of erythropoiesis, as consideration of healthy person's state: 
- Initial erythropoiesis species composition
	
	$E_{init} = \frac{Hkt_{init} \cdot Vol_{blood}}{Vol_{E} + \frac{t_{R,a,init}}{t_{E,d}} \cdot Vol_{R}}$
	$R_{init} = E_{init} \cdot \frac{t_{R,a,init}}{t_{E,d}}$, here lies one problem why RPI low $R_{init}<0.01\cdot E_{init}$
	$P_{init} = \frac{1}{1024}\cdot R_{init}\cdot \frac{t_{P,a,init}} {t_{R,a,init}}$
	$LDH_{Blood}=\frac{[LDH] \cdot \frac{ln(2)} {t_{hlf,LDH,d}} \cdot Vol_{blood}}{[E]\cdot \frac{2 \cdot ln(2)}{t_{E,d}} +\underbrace{J_{R,d}}_{=0}}$ 

- Birth rate of precursor cells, calculated from steady state and then set for whole time 
	$k_{P,b}= P\cdot \frac{( 1 + k_{P,d} * Hb^{r_{P,d}})} {a_{P,d}}+P(\frac{2\cdot ln(2)}{t_{P,a}})$ 
	$Hb_{init} = \frac{(Vol_E \cdot E_{init} * Hb_{conc,E} + Vol_R \cdot R_{init} \cdot Hb_{conc,R})}{(10 \cdot Vol_{blood})}$ in model change of concentration unit from g/l to g/dl, thus div. by 10