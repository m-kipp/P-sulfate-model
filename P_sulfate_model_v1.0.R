##### P-SO4 model
##### Michael Kipp
##### May 2022 

## modern flux & reservoir values - read in these values before executing code for any of the figures below

P_modern=2900e12 # moles of P in modern ocean; from Ruttenberg (2003) Table 1
O2_modern=38e18 # moles O2 in modern atmosphere
SO4_modern=40e18 # moles SO4 in modern ocean
O2_weathering_modern=3.335761e13 # modern oxidative weathering flux, mol/yr
NPP_modern=4e15 # modern net primary productivity; mol/yr, from Field et al (1999)
SO4_input_modern=4.5e12 # modern riverine flux of SO4 to ocean; mol/yr
P_input_modern=0.27e12 # modern riverine flux of P to ocean; mol/yr of "reactive" (i.e. bioavailable, not particulate) from Ruttenberg (2003) Table 3
Redfield_CP=106 # molar C/P ratio of organic matter; Redfield (1958)
SO4_hsc=0.00333 # half-saturation constant for SO4 reduction (0.00333=100uM, as in Pallud & Van Cappellen 2006)
O2_hsc=1e-4 # from Slomp and Van Cappellen (2007) Biogeosciences; Alcott et al (2019) Science
e_b=0.00722 # calibrated for steady state; in general agreement with Holland (1984); Hedges & Keil (1995); Middelburg (2019)
f_sulfate_modern=0.5 # from S isotope record, e.g. Canfield + Farquhar (2009) PNAS
gypsum_hsc=0.2 # tuned to give negligible gypsum burial in Archean
Anoxic_CP_modern=3650 # 4000 in Van Cappellen & Ingall (1996) and Lenton et al (2017); 1100 in Alcott et al (2019) and Slomp & Van Cappellen (2007)
reduced_outgassing=0 # no additional O2 sink in base model 

### Model run at modern steady state values for 1 Myr. 

i=2
max_time=1e6 # length of run; yr
time_step=1e3 # time step; yr
t=seq(0,max_time,by=time_step)
P=NULL
P[1]=P_modern # initial marine P inventory

O2=NULL
O2[1]=O2_modern # initial atmospheric O2 inventory 
SO4=NULL
SO4[1]=SO4_modern # initial marine SO4 inventory

out=NULL

for (i in 2:length(t)){
	
	time=t[i]
	
	P_input=P_input_modern*time_step
	
	Corg_prod=(P[i-1]/P_modern)*NPP_modern*time_step 
	
	resp_O2=((O2[i-1]/((O2_hsc*O2_modern)+O2[i-1])))
	resp_SO4=((SO4[i-1]/((SO4_hsc*SO4_modern)+SO4[i-1])))
	total_pot_resp=resp_O2+resp_SO4
	if(total_pot_resp>0.99278){total_e_b=1-0.99278} else{total_e_b=1-total_pot_resp}
	
	Corg_burial= total_e_b* Corg_prod
	
	f_anoxic=1/(1+exp(-12*(0.5*(P[i-1]/P_modern)-(O2[i-1]/O2_modern)))) # eqn & coefficients from Lenton et al (2018) Earth Sci Rev

	SO4_input=sqrt(O2[i-1]/O2_modern)*SO4_input_modern*time_step
	
	J_sulfate=f_sulfate_modern*((SO4[i-1]/((gypsum_hsc*SO4_modern)+SO4[i-1]))*SO4_input_modern*time_step)
	J_sulfide=(1-f_sulfate_modern)*((SO4[i-1]/((SO4_hsc*SO4_modern)+SO4[i-1]))*SO4_input_modern*time_step)
	SO4[i]=SO4[i-1]+SO4_input-J_sulfide-J_sulfate
	
	P_oxic_burial=(1-f_anoxic)*Corg_burial/Redfield_CP
	Anoxic_CP=if((Anoxic_CP_modern*(SO4[i-1]/((SO4_hsc*SO4_modern)+SO4[i-1])))>106){(Anoxic_CP_modern*(SO4[i-1]/((SO4_hsc*SO4_modern)+SO4[i-1])))} else{106}
	P_anoxic_burial=(f_anoxic)*Corg_burial/Anoxic_CP
	P_total_burial= P_oxic_burial+ P_anoxic_burial
	
	P[i]=if(P[i-1]+P_input-P_total_burial>0){P[i-1]+P_input-P_total_burial} else{0}
	
	O2_source=Corg_burial+2*J_sulfide
	O2_sink=(O2_weathering_modern*time_step)*sqrt(O2[i-1]/O2_modern)+reduced_outgassing*time_step
	O2[i]=if(O2[i-1]+O2_source-O2_sink>0){O2[i-1]+O2_source-O2_sink} else{0}
		
	results=c(time,P[i], P_input/time_step, P_oxic_burial/time_step, P_anoxic_burial/time_step, P_total_burial/time_step, f_anoxic,O2[i], P_oxic_burial/P_total_burial,P_anoxic_burial/P_total_burial, Corg_burial/time_step, Anoxic_CP,SO4[i], Corg_burial/Corg_prod, SO4_input/time_step,J_sulfate/(J_sulfate+J_sulfide), O2_source/time_step,O2_sink/time_step)
	out=rbind(out,results)
	
}

colnames(out)=c("time","P_ocean","P_input","P_oxic_burial","P_anoxic_burial","P_total_burial","f_anoxic","O2_atmosphere","rel_oxic_P_burial","rel_anoxic_P_burial","Corg_burial","CP_anoxic","SO4_ocean","e_b","SO4_input","f_sulfate","O2_source","O2_sink")

par(mfrow=c(3,3))
par(mar=c(4,4,1,1))
par(oma=c(0,0,0,0))
plot(out[,2]~out[,1],cex=0,ylim=c(2000e9,6000e12),xlab=NA,ylab=NA,axes=F,log="y")
lines(out[,2]~out[,1])
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(2.9e15,2.9e14,2.9e13,2.9e12),labels=c(100,10,1,0.1))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"[P] (% modern)")

plot(out[,6]~out[,1],type="l",xlab=NA,ylab=NA,axes=F,ylim=c(2e10,5e14),log="y")
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(1e10,1e11,1e12,1e13,1e14,1e15),labels=c("1e10","1e11","1e12","1e13","1e14","1e15"))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"total P burial (mol/yr)")

plot(out[,7]~out[,1],type="l",xlab=NA,ylab=NA,axes=F,ylim=c(0.001,1),log="y")
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(0.001,0.01,0.1,1),labels=c("0.001","0.01","0.1","1"))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"f_anoxic (%)")

plot(out[,8]~out[,1],type="l",ylim=c(30e12,80e18),xlab=NA,ylab=NA,axes=F,log="y")
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(38e18,1e-2*38e18,1e-4*38e18,1e-6*38e18),labels=c("1e0","1e-2","1e-4","1e-6"))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"pO2 (PAL)")

plot(out[,11]~out[,1],type="l",xlab=NA,ylab=NA,axes=F,ylim=c(2e12,5e16),log="y")
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(1e13,1e14,1e15,1e16),labels=c("1e13","1e14","1e15","1e16"))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"C_org_burial (mol/yr)")

plot(out[,13]~out[,1],type="l",xlab=NA,ylab=NA,axes=F,ylim=c(1e15,1e20),log="y")
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(40e18,40e17,40e16,40e15,40e14),labels=c("30mM","3mM","300uM","30uM","3uM"))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"[SO4]")

plot(out[,15]~out[,1],type="l",xlab=NA,ylab=NA,axes=F,ylim=c(1e11,1e14),log="y")
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(1e11,1e12,1e13,1e14,1e15,1e16,1e17,1e18),labels=c("1e11","1e12","1e13","1e14","1e15","1e16","1e17","1e18"))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"SO4 input (mol/yr)")

plot(out[,14]~out[,1],type="l",xlab=NA,ylab=NA,axes=F,ylim=c(0.001,1),log="y")
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(0.001,0.01,0.1,1),labels=c("0.001","0.01","0.1","1"))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"e_b")

plot(out[,16]~out[,1],type="l",xlab=NA,ylab=NA,axes=F,ylim=c(0,1))
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"f_sulfate")







### Figure 2. Long-term negative feedback on pO2. 

i=2
max_time=20e6 # length of run; yr
time_step=1e4 # time step; yr
t=seq(0,max_time,by=time_step)
P=NULL
P[1]=P_modern # initial marine P inventory

O2=NULL
O2[1]=O2_modern # initial atmospheric O2 inventory 
SO4=NULL
SO4[1]=SO4_modern*1e-4 # initial marine SO4 inventory; *1e-4 for low-SO4 test
O2_weathering_modern=3.335761e13*(0.865576) # *0.865576 for low-SO4 test (since sulfide O2 source becomes negligible)

out=NULL

reduced_outgassing_seq=seq(0,0,length.out=length(t))
reduced_outgassing_seq[200:length(t)]= O2_weathering_modern*0.25 # 25% increase in O2 sink

for (i in 2:length(t)){
	
	time=t[i]
	
	P_input=P_input_modern*time_step
	
	Corg_prod=(P[i-1]/P_modern)*NPP_modern*time_step 
	
	resp_O2=((O2[i-1]/((O2_hsc*O2_modern)+O2[i-1])))
	resp_SO4=((SO4[i-1]/((SO4_hsc*SO4_modern)+SO4[i-1])))
	total_pot_resp=resp_O2+resp_SO4
	if(total_pot_resp>0.99278){total_e_b=1-0.99278} else{total_e_b=1-total_pot_resp}
	
	Corg_burial= total_e_b* Corg_prod
	
	f_anoxic=1/(1+exp(-12*(0.5*(P[i-1]/P_modern)-(O2[i-1]/O2_modern)))) # eqn & coefficients from Lenton et al (2018) Earth Sci Rev

	SO4_input=sqrt(O2[i-1]/O2_modern)*SO4_input_modern*time_step*1e-2 # *1e-2 for low-SO4 test
	
	J_sulfate=f_sulfate_modern*((SO4[i-1]/((gypsum_hsc*SO4_modern)+SO4[i-1]))*SO4_input_modern*time_step)
	J_sulfide=(1-f_sulfate_modern)*((SO4[i-1]/((SO4_hsc*SO4_modern)+SO4[i-1]))*SO4_input_modern*time_step)
	SO4[i]=SO4[i-1]+SO4_input-J_sulfide-J_sulfate
	
	P_oxic_burial=(1-f_anoxic)*Corg_burial/Redfield_CP
	Anoxic_CP=if((Anoxic_CP_modern*(SO4[i-1]/((SO4_hsc*SO4_modern)+SO4[i-1])))>106){(Anoxic_CP_modern*(SO4[i-1]/((SO4_hsc*SO4_modern)+SO4[i-1])))} else{106}
	P_anoxic_burial=(f_anoxic)*Corg_burial/Anoxic_CP
	P_total_burial= P_oxic_burial+ P_anoxic_burial
	
	P[i]=if(P[i-1]+P_input-P_total_burial>0){P[i-1]+P_input-P_total_burial} else{0}
	
	O2_source=Corg_burial+2*J_sulfide
	O2_sink=(O2_weathering_modern*time_step)*sqrt(O2[i-1]/O2_modern)+ reduced_outgassing_seq[i]*time_step
	O2[i]=if(O2[i-1]+O2_source-O2_sink>0){O2[i-1]+O2_source-O2_sink} else{0}
		
	results=c(time,P[i], P_input/time_step, P_oxic_burial/time_step, P_anoxic_burial/time_step, P_total_burial/time_step, f_anoxic,O2[i], P_oxic_burial/P_total_burial,P_anoxic_burial/P_total_burial, Corg_burial/time_step, Anoxic_CP,SO4[i], Corg_burial/Corg_prod, SO4_input/time_step,J_sulfate/(J_sulfate+J_sulfide), O2_source/time_step,O2_sink/time_step)
	out=rbind(out,results)
	
}

colnames(out)=c("time","P_ocean","P_input","P_oxic_burial","P_anoxic_burial","P_total_burial","f_anoxic","O2_atmosphere","rel_oxic_P_burial","rel_anoxic_P_burial","Corg_burial","CP_anoxic","SO4_ocean","e_b","SO4_input","f_sulfate","O2_source","O2_sink")

par(mfrow=c(3,4))
par(mar=c(4,4,1,1))
par(oma=c(0,0,0,0))
plot(out[,2]~out[,1],cex=0,ylim=c(2000e9,6000e12),xlab=NA,ylab=NA,axes=F,log="y")
lines(out[,2]~out[,1])
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(2.9e15,2.9e14,2.9e13,2.9e12),labels=c(100,10,1,0.1))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"[P] (% modern)")

plot(out[,6]~out[,1],type="l",xlab=NA,ylab=NA,axes=F,ylim=c(2e10,5e14),log="y")
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(1e10,1e11,1e12,1e13,1e14,1e15),labels=c("1e10","1e11","1e12","1e13","1e14","1e15"))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"total P burial (mol/yr)")

plot(out[,7]~out[,1],type="l",xlab=NA,ylab=NA,axes=F,ylim=c(0.001,1),log="y")
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(0.001,0.01,0.1,1),labels=c("0.001","0.01","0.1","1"))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"f_anoxic (%)")

plot(out[,8]~out[,1],type="l",ylim=c(30e12,80e18),xlab=NA,ylab=NA,axes=F)
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(38e18,1e-2*38e18,1e-4*38e18,1e-6*38e18),labels=c("1e0","1e-2","1e-4","1e-6"))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"pO2 (PAL)")

plot(out[,11]~out[,1],type="l",xlab=NA,ylab=NA,axes=F,ylim=c(2e12,5e16),log="y")
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(1e13,1e14,1e15,1e16),labels=c("1e13","1e14","1e15","1e16"))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"C_org_burial (mol/yr)")

plot(out[,13]~out[,1],type="l",xlab=NA,ylab=NA,axes=F,ylim=c(1e15,1e20),log="y")
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(40e18,40e17,40e16,40e15,40e14),labels=c("30mM","3mM","300uM","30uM","3uM"))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"[SO4]")

plot(out[,15]~out[,1],type="l",xlab=NA,ylab=NA,axes=F,ylim=c(1e11,1e14),log="y")
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(1e11,1e12,1e13,1e14,1e15,1e16,1e17,1e18),labels=c("1e11","1e12","1e13","1e14","1e15","1e16","1e17","1e18"))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"SO4 input (mol/yr)")

plot(out[,14]~out[,1],type="l",xlab=NA,ylab=NA,axes=F,ylim=c(0.001,1),log="y")
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(0.001,0.01,0.1,1),labels=c("0.001","0.01","0.1","1"))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"e_b")

plot(out[,16]~out[,1],type="l",xlab=NA,ylab=NA,axes=F,ylim=c(0,1))
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"f_sulfate")

plot(out[,17]~out[,1],type="l",xlab=NA,ylab=NA,axes=F,ylim=c(1e13,1e14),log="y")
lines(out[,18]~out[,1],col="red",lty=2)
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"O2_flux")


### Figure 3. OAE amplification.  

i=2
max_time=1e6 # length of run; yr
time_step=1e3 # time step; yr
t=seq(0,max_time,by=time_step)
P=NULL
P[1]=P_modern # initial marine P inventory

O2=NULL
O2[1]=O2_modern # initial atmospheric O2 inventory 
SO4=NULL
SO4[1]=SO4_modern*1e-4 # initial marine SO4 inventory; *1e-4 for low-SO4 test

SO4_input_modern=4.5e12
O2_weathering_modern=3.335761e13*(0.865576) # *0.865576 for low-SO4 test (since sulfide O2 source becomes negligible)
reduced_outgassing=0 

out=NULL

P_input_seq=seq(P_input_modern,P_input_modern,length.out=length(t))
P_input_seq[200:300]= P_input_modern*1.5 # 50% increase in P input

for (i in 2:length(t)){
	
	time=t[i]
	
	P_input= P_input_seq[i]*time_step
	
	Corg_prod=(P[i-1]/P_modern)*NPP_modern*time_step 
	
	resp_O2=((O2[i-1]/((O2_hsc*O2_modern)+O2[i-1])))
	resp_SO4=((SO4[i-1]/((SO4_hsc*SO4_modern)+SO4[i-1])))
	total_pot_resp=resp_O2+resp_SO4
	if(total_pot_resp>0.99278){total_e_b=1-0.99278} else{total_e_b=1-total_pot_resp}
	
	Corg_burial= total_e_b* Corg_prod
	
	f_anoxic=1/(1+exp(-12*(0.5*(P[i-1]/P_modern)-(O2[i-1]/O2_modern)))) # eqn & coefficients from Lenton et al (2018) Earth Sci Rev

	SO4_input=sqrt(O2[i-1]/O2_modern)*SO4_input_modern*time_step*1e-2 # *1e-2 for low-SO4 test
	
	J_sulfate=f_sulfate_modern*((SO4[i-1]/((gypsum_hsc*SO4_modern)+SO4[i-1]))*SO4_input_modern*time_step)
	J_sulfide=(1-f_sulfate_modern)*((SO4[i-1]/((SO4_hsc*SO4_modern)+SO4[i-1]))*SO4_input_modern*time_step)
	SO4[i]=SO4[i-1]+SO4_input-J_sulfide-J_sulfate
	
	P_oxic_burial=(1-f_anoxic)*Corg_burial/Redfield_CP
	Anoxic_CP=if((Anoxic_CP_modern*(SO4[i-1]/((SO4_hsc*SO4_modern)+SO4[i-1])))>106){(Anoxic_CP_modern*(SO4[i-1]/((SO4_hsc*SO4_modern)+SO4[i-1])))} else{106}
	P_anoxic_burial=(f_anoxic)*Corg_burial/Anoxic_CP
	P_total_burial= P_oxic_burial+ P_anoxic_burial
	
	P[i]=if(P[i-1]+P_input-P_total_burial>0){P[i-1]+P_input-P_total_burial} else{0}
	
	O2_source=Corg_burial+2*J_sulfide
	O2_sink=(O2_weathering_modern*time_step)*sqrt(O2[i-1]/O2_modern)+ reduced_outgassing*time_step
	O2[i]=if(O2[i-1]+O2_source-O2_sink>0){O2[i-1]+O2_source-O2_sink} else{0}
		
	results=c(time,P[i], P_input/time_step, P_oxic_burial/time_step, P_anoxic_burial/time_step, P_total_burial/time_step, f_anoxic,O2[i], P_oxic_burial/P_total_burial,P_anoxic_burial/P_total_burial, Corg_burial/time_step, Anoxic_CP,SO4[i], Corg_burial/Corg_prod, SO4_input/time_step,J_sulfate/(J_sulfate+J_sulfide), O2_source/time_step,O2_sink/time_step)
	out=rbind(out,results)
	
}

colnames(out)=c("time","P_ocean","P_input","P_oxic_burial","P_anoxic_burial","P_total_burial","f_anoxic","O2_atmosphere","rel_oxic_P_burial","rel_anoxic_P_burial","Corg_burial","CP_anoxic","SO4_ocean","e_b","SO4_input","f_sulfate","O2_source","O2_sink")

par(mfrow=c(3,4))
par(mar=c(4,4,1,1))
par(oma=c(0,0,0,0))
plot(out[,2]~out[,1],cex=0,ylim=c(2000e9,6000e12),xlab=NA,ylab=NA,axes=F,log="y")
lines(out[,2]~out[,1])
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(2.9e15,2.9e14,2.9e13,2.9e12),labels=c(100,10,1,0.1))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"[P] (% modern)")

plot(out[,6]~out[,1],type="l",xlab=NA,ylab=NA,axes=F,ylim=c(2e10,5e14),log="y")
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(1e10,1e11,1e12,1e13,1e14,1e15),labels=c("1e10","1e11","1e12","1e13","1e14","1e15"))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"total P burial (mol/yr)")

plot(out[,7]~out[,1],type="l",xlab=NA,ylab=NA,axes=F,ylim=c(0.001,1),log="y")
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(0.001,0.01,0.1,1),labels=c("0.001","0.01","0.1","1"))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"f_anoxic (%)")

plot(out[,8]~out[,1],type="l",ylim=c(30e12,80e18),xlab=NA,ylab=NA,axes=F)
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(38e18,1e-2*38e18,1e-4*38e18,1e-6*38e18),labels=c("1e0","1e-2","1e-4","1e-6"))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"pO2 (PAL)")

plot(out[,11]~out[,1],type="l",xlab=NA,ylab=NA,axes=F,ylim=c(2e12,5e16),log="y")
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(1e13,1e14,1e15,1e16),labels=c("1e13","1e14","1e15","1e16"))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"C_org_burial (mol/yr)")

plot(out[,13]~out[,1],type="l",xlab=NA,ylab=NA,axes=F,ylim=c(1e15,1e20),log="y")
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(40e18,40e17,40e16,40e15,40e14),labels=c("30mM","3mM","300uM","30uM","3uM"))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"[SO4]")

plot(out[,15]~out[,1],type="l",xlab=NA,ylab=NA,axes=F,ylim=c(1e11,1e14),log="y")
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(1e11,1e12,1e13,1e14,1e15,1e16,1e17,1e18),labels=c("1e11","1e12","1e13","1e14","1e15","1e16","1e17","1e18"))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"SO4 input (mol/yr)")

plot(out[,14]~out[,1],type="l",xlab=NA,ylab=NA,axes=F,ylim=c(0.001,1),log="y")
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(0.001,0.01,0.1,1),labels=c("0.001","0.01","0.1","1"))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"e_b")

plot(out[,16]~out[,1],type="l",xlab=NA,ylab=NA,axes=F,ylim=c(0,1))
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"f_sulfate")

plot(out[,17]~out[,1],type="l",xlab=NA,ylab=NA,axes=F,ylim=c(1e13,1e14),log="y")
lines(out[,18]~out[,1],col="red",lty=2)
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"O2_flux")



### Figure 4. Whiff of O2.  

i=2
max_time=600e3 # length of run; yr
time_step=50 # time step; yr
t=seq(0,max_time,by=time_step)
P=NULL
P[1]=P_modern*1e-2 # initial marine P inventory; *1e-2 for Archean P reservoir

O2=NULL
O2[1]=O2_modern*1e-6 # initial atmospheric O2 inventory; *1e-6 for Archean pO2
SO4=NULL
SO4[1]=SO4_modern*1e-4 # initial marine SO4 inventory; *1e-4 for Archean SO4 reservoir

SO4_input_volcanic=0.5e11 # Archean SO4 source via photolysis of volcanic SO2
SO4_input_volcanic_seq=seq(SO4_input_volcanic, SO4_input_volcanic,length.out=length(t))
#SO4_input_volcanic_seq[2000:3000]= SO4_input_volcanic*5 # alt. test: 5x increase in volcanic SO2

P_input_seq=seq(P_input_modern,P_input_modern,length.out=length(t))
P_input_seq[3000:4000]= P_input_modern*1.5 # 40-50% increase in P input

out=NULL

reduced_outgassing=40e12 # Kipp et al (2021): 13.5e12 baseline, 24.3e12 (0.5 log unit fO2 change), 43.1e12 (1.0 log unit fO2 change)

for (i in 2:length(t)){
	
	time=t[i]
	
	P_input= P_input_seq[i]*time_step
	
	Corg_prod=(P[i-1]/P_modern)*NPP_modern*time_step 
	
	resp_O2=((O2[i-1]/((O2_hsc*O2_modern)+O2[i-1])))
	resp_SO4=((SO4[i-1]/((SO4_hsc*SO4_modern)+SO4[i-1])))
	total_pot_resp=resp_O2+resp_SO4
	if(total_pot_resp>0.99278){total_e_b=1-0.99278} else{total_e_b=1-total_pot_resp}
	
	Corg_burial= total_e_b* Corg_prod
	
	f_anoxic=1/(1+exp(-12*(0.5*(P[i-1]/P_modern)-(O2[i-1]/O2_modern)))) # eqn & coefficients from Lenton et al (2018) Earth Sci Rev

	SO4_input=sqrt(O2[i-1]/O2_modern)*SO4_input_modern*time_step+ SO4_input_volcanic_seq[i]*time_step
	
	J_sulfate=f_sulfate_modern*((SO4[i-1]/((gypsum_hsc*SO4_modern)+SO4[i-1]))*SO4_input_modern*time_step)
	J_sulfide=(1-f_sulfate_modern)*((SO4[i-1]/((SO4_hsc*SO4_modern)+SO4[i-1]))*SO4_input_modern*time_step)
	SO4[i]=SO4[i-1]+SO4_input-J_sulfide-J_sulfate
	
	P_oxic_burial=(1-f_anoxic)*Corg_burial/Redfield_CP
	Anoxic_CP=106 #if((Anoxic_CP_modern*(SO4[i-1]/((SO4_hsc*SO4_modern)+SO4[i-1])))>106){(Anoxic_CP_modern*(SO4[i-1]/((SO4_hsc*SO4_modern)+SO4[i-1])))} else{106} # =106 for "control" run with no anoxic P recycling
	P_anoxic_burial=(f_anoxic)*Corg_burial/Anoxic_CP
	P_total_burial= P_oxic_burial+ P_anoxic_burial
	Total_CP=(P_anoxic_burial/P_total_burial)*Anoxic_CP+(P_oxic_burial/P_total_burial)*Redfield_CP
	
	P[i]=if(P[i-1]+P_input-P_total_burial>0){P[i-1]+P_input-P_total_burial} else{0}
	
	O2_source=Corg_burial+2*J_sulfide
	O2_sink=(O2_weathering_modern*time_step)*sqrt(O2[i-1]/O2_modern)+ reduced_outgassing*time_step
	O2[i]=if(O2[i-1]+O2_source-O2_sink>0){O2[i-1]+O2_source-O2_sink} else{0}
		
	results=c(time,P[i], P_input/time_step, P_oxic_burial/time_step, P_anoxic_burial/time_step, P_total_burial/time_step, f_anoxic,O2[i], P_oxic_burial/P_total_burial,P_anoxic_burial/P_total_burial, Corg_burial/time_step, Anoxic_CP,SO4[i], Corg_burial/Corg_prod, SO4_input/time_step,J_sulfate/(J_sulfate+J_sulfide), O2_source/time_step,O2_sink/time_step,Total_CP)
	out=rbind(out,results)
	
}

colnames(out)=c("time","P_ocean","P_input","P_oxic_burial","P_anoxic_burial","P_total_burial","f_anoxic","O2_atmosphere","rel_oxic_P_burial","rel_anoxic_P_burial","Corg_burial","CP_anoxic","SO4_ocean","e_b","SO4_input","f_sulfate","O2_source","O2_sink","CP_total")


par(mfrow=c(3,4))
par(mar=c(4,4,1,1))
par(oma=c(0,0,0,0))
plot(out[,2]~out[,1],cex=0,ylim=c(2000e9,6000e12),xlab=NA,ylab=NA,axes=F,log="y")
lines(out[,2]~out[,1])
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(2.9e15,2.9e14,2.9e13,2.9e12),labels=c(100,10,1,0.1))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"[P] (% modern)")

plot(out[,6]~out[,1],type="l",xlab=NA,ylab=NA,axes=F,ylim=c(2e10,5e14),log="y")
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(1e10,1e11,1e12,1e13,1e14,1e15),labels=c("1e10","1e11","1e12","1e13","1e14","1e15"))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"total P burial (mol/yr)")

plot(out[,7]~out[,1],type="l",xlab=NA,ylab=NA,axes=F,ylim=c(0.001,1),log="y")
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(0.001,0.01,0.1,1),labels=c("0.001","0.01","0.1","1"))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"f_anoxic (%)")

plot(out[,8]~out[,1],type="l",ylim=c(30e12,80e18),xlab=NA,ylab=NA,axes=F,log="y")
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(38e18,1e-2*38e18,1e-4*38e18,1e-6*38e18),labels=c("1e0","1e-2","1e-4","1e-6"))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"pO2 (PAL)")

plot(out[,11]~out[,1],type="l",xlab=NA,ylab=NA,axes=F,ylim=c(2e12,5e16),log="y")
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(1e13,1e14,1e15,1e16),labels=c("1e13","1e14","1e15","1e16"))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"C_org_burial (mol/yr)")

plot(out[,13]~out[,1],type="l",xlab=NA,ylab=NA,axes=F,ylim=c(1e15,1e20),log="y")
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(40e18,40e17,40e16,40e15,40e14),labels=c("30mM","3mM","300uM","30uM","3uM"))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"[SO4]")

plot(out[,15]~out[,1],type="l",xlab=NA,ylab=NA,axes=F,ylim=c(1e10,1e14),log="y")
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(1e11,1e12,1e13,1e14,1e15,1e16,1e17,1e18),labels=c("1e11","1e12","1e13","1e14","1e15","1e16","1e17","1e18"))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"SO4 input (mol/yr)")

plot(out[,14]~out[,1],type="l",xlab=NA,ylab=NA,axes=F,ylim=c(0.001,1),log="y")
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(0.001,0.01,0.1,1),labels=c("0.001","0.01","0.1","1"))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"e_b")

plot(out[,16]~out[,1],type="l",xlab=NA,ylab=NA,axes=F,ylim=c(0,1))
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"f_sulfate")

plot(out[,17]~out[,1],type="l",xlab=NA,ylab=NA,axes=F,ylim=c(1e13,1e14),log="y")
lines(out[,18]~out[,1],col="red",lty=2)
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"O2_flux")

plot(out[,19]~out[,1],type="l",xlab=NA,ylab=NA,axes=F,ylim=c(50,200))
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"CP_total")

plot(out[,3]~out[,1],type="l",xlab=NA,ylab=NA,axes=F,ylim=c(2e10,5e12),log="y")
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(1e10,1e11,1e12,1e13,1e14,1e15),labels=c("1e10","1e11","1e12","1e13","1e14","1e15"))
box(lwd=1)
mtext(side=1,line=2,"time (yr)")
mtext(side=2,line=2.5,"P input (mol/yr)")


