

#R3.6.2
require(MASS);

Trun=1500; nReps=100;
#Baseline parameter sets following scenrios in panels a-d of Figure 2.
#Some parameter name differences here compared to manuscript: Acor=R, rLt=r_T, rPt=p_T, rP=p_M, rPc=p_NM, F=f, xi=kappa, Glo=E
Pmat=rbind( c(L=8,N=18,c=0.8,Alph=0,Acor=0.0,rL=0.0,rLt=0.4,rPt=0.5,rP=4,rPc=1.6,F=0,mu=5/1e2,beta=.75,xi=0.25,Glo=1,f0=0.15,T=0.5),
			c(L=8,N=18,c=0.8,Alph=0,Acor=0.0,rL=0.6,rLt=0.4,rPt=0.5,rP=4,rPc=1.6,F=0,mu=5/1e2,beta=.75,xi=0.25,Glo=1,f0=0.15,T=0.5),
			c(L=8,N=18,c=0.8,Alph=2,Acor=0.7,rL=0.0,rLt=0.4,rPt=0.5,rP=4,rPc=1.6,F=0,mu=5/1e2,beta=.75,xi=0.25,Glo=1,f0=0.15,T=0.5),
			c(L=8,N=18,c=0.8,Alph=2,Acor=0.7,rL=0.6,rLt=0.4,rPt=0.5,rP=4,rPc=1.6,F=0,mu=5/1e2,beta=.75,xi=0.25,Glo=1,f0=0.15,T=0.5))

#Below are parameter sets for various simulations. A loop at the end runs simulations over chosen set of conditions Starts			
			
#### Simlations for Figures 1 and 2:
Starts=as.matrix(expand.grid(Pset=1:4,c=0.8,T=c(0,0.5,1.5),cvar=0,Gcor=0,Glap=0,infl=0,Rep=1:nReps))

#### Simlations for Figure 3 - Effect of constitutive punishment costs on overlapping institutions
# Note that agreement overlap (parameter b) here is implemented as 2 parameters:
# T remains the global mitigation frequency m_G at which global agreement attempts begin,
# and Tlocoff is now the m_g level at which local agreements are phased out
Fs=seq(0,0.65,len=30); 
Starts=rbind(as.matrix(expand.grid(Pset=4,c=0.8,T=c(0,0.5),cvar=0,Gcor=0,Glap=0,infl=0,F=Fs,Tlocoff=c(0.5,1.5),Rep=1:nReps)),
	as.matrix(expand.grid(Pset=4,c=0.8,T=c(0,0.5),cvar=0,Gcor=0,Glap=0,infl=0,F=Fs,Tlocoff=c(0,1.5),Rep=1:nReps)));
Starts=Starts[-which(Starts[,"T"]==1.5 & Starts[,"Tlocoff"]==0),]


#### Simulations for Figure 4 - effect of inequality in emissions and negotiations
# Some parameter name differences here compared to Methods in manuscript: infl=rho, Gcor=rho_G
# Glap denotes overlap in group/local agreement membership, with Glap=0 used in the main analysis (ie no group overlap)
# and Glap=1 corresponding to every player belonging to every local group (equivalent to a single global agreement)
ntests=10; lapSeq=c(0,round(exp(-seq(4.5,2,len=ntests-1)),3)); inflSeq=seq(0.01,1,len=ntests);
Starts=rbind(expand.grid(Pset=4,c=0.8,T=c(0,0.5,1.5),cvar=c(0,6),Gcor=c(0,0.85),Glap=0,infl=inflSeq,Rep=1:nReps),
	expand.grid(Pset=4,c=0.8,T=c(0,0.5,1.5),cvar=c(0,6),Gcor=seq(0,0.85,len=ntests),Glap=0,infl=c(0.1,0.2,0.4,0.6),Rep=1:nReps),
	expand.grid(Pset=4,c=0.8,T=c(0,0.5,1.5),cvar=c(0,6),Gcor=0,Glap=lapSeq,infl=c(0.1,0.2,0.4,0.6),Rep=1:nReps))







#Agreement switching function Theta
sig=function(x,thresh,smoo=15) 1-1/(1+exp(smoo*(x-thresh)))

#Function to generate lognormal distribution from set of quantiles
qlnormT=function(q,mu,sig) qlnorm(q, meanlog=log(mu/sqrt(1+(sig/mu)^2)), sdlog=sqrt(log(1+(sig/mu)^2)) )

#Function assigning players to groups and the emissions of each player
group=function(parms,seed0){
	ss0=function(S=seed0) set.seed(S); L=parms["L"]; np=parms["N"]*L; G=matrix(0,np,L); colnames(G)=paste0("G",1:L);
	ss0(); qts=mvrnorm(n=np,mu=c(0,0), diag(2)*(1-parms["Gcor"])+parms["Gcor"]);
	Cs=round(qlnormT(pnorm(qts[,2]),parms["c"],parms["cvar"]),3); ps=0.5+parms["Gcor"]*(0.5 - (1:L)/L);
	ss0(); Ns=rmultinom(20,np,ps); Ns=Ns[,c(which(Ns[L,]>1),1)[1]]; #cor(Cs,rep(Ns,Ns)[rank(qts[,1])]); #cor(1:L,ps)=1 => final cor := Gcor in qts
	G[cbind(1:np, rep(1:L,Ns)[rank(qts[,1])])]=1; ss0(); G[G==0]=rbinom(np*(L-1),1,parms["Glap"]); #Glap=degree of grp overlap (=prop multi-group memberships realized)
	ss0(); cbind(Cs,M=Cs*rbinom(np,1,parms["f0"]),G);
}

#Function calculating payoffs for each step
payfun3=function(pS,cS,parms,pSG,N=parms["N"],Li=length(pS)){
	#pS is fraction of players mitigating; cS is either cost for deciding player OR (for group rivalry) mean cost of currently mitigating players in each group
	#Determine at what scales punishment is happening, ultimately calculating punishment fines PF
	pSL=pS; if(parms["infl"]>0){ pSL=cS; cS=parms["c"]; };
	Ps=sig(c(pSG,pSG,pSG,pSL),rep(parms[c("T","Tlocoff","rPt")],c(1,1,1+Li)));
	Attempt=(1-Ps[2])+Ps[1]; dlt=(1-Ps[2])*Ps[-(1:3)] + parms["Glo"]*Ps[1]*Ps[3];
	#Effects of punishment PF on each strategy in each group
	PF=-cbind(parms["F"]*Attempt + (1-parms["F"])*dlt, dlt)%*%diag(parms[c("rPc","rP")])
	#Add mitigator costs - which can vary among groups due to economies of scale
	Pays = PF - cbind(cS*(1-parms["rL"]*sig(pS,parms["rLt"])),  0)
	if(Li==1) return(as.vector(Pays)); return(Pays[,1]*pS + Pays[,2]*(1-pS)); 
}

#Function to decide whether a chosen player switches their climate policy
m9iter=function(S,dSi,parms,seed){ #Threshold punishment model with only 2 strategies and chickens
	#Start off by calculating group-level stats, seperately tracking pS and mean cost of mitigators.
	#belonging denotes the groups to which the player reconsidering policy belongs (single integer for Glap=0)
	L=parms["L"]; infl=parms["infl"]; Si=S[dSi,]; belonging=(1:L)[Si[-(1:2)]==1]; Mitigator=Si[2]>0;
	v=1+0*S[,1]; SS=cbind(Nj=(v%*%S[,-(1:2)])[1,], nM=(v%*%((S[,2]>0)*S[,-(1:2)]))[1,], sumC=(v%*%(S[,2]*S[,-(1:2)]))[1,]); 
	SS=cbind(SS, pS=SS[,"nM"]/SS[,"Nj"], avC=SS[,"sumC"]/SS[,"nM"]); SS[is.na(SS)]=0; #NA mean costs are actually zero

	Pays=matrix(nrow=0,ncol=2); pSG=mean(S[,2]>0); cS1s=rep(Si[1],L);
	if(infl>0){ SGC=(v%*%(S[,1]*S[,-(1:2)]))[1,]; grpS=cS1s=(1-infl)*SS[,"pS"] + infl*SS[,"sumC"]/SGC; }
	SGA=(1-infl)*(mean(S[,2]>0)+(1:-1)/nrow(S)) + infl*(sum(S[,2])+(1:-1)*Si[1])/sum(S[,1])
	for(i in belonging) Pays=rbind(Pays,payfun3(SS[i,2]/SS[i,1],cS1s[i],parms,SGA[2])); Pays=colMeans(Pays);
	dPg=0; if(parms["Alph"]>0){ dPg=1[0]; for(k in belonging){
			pSk=matrix(SS[,"pS"],L,3); pSk[k,]=pSk[k,]+(1:-1)/SS[k,1]; 
			if(infl==0){ cSk=matrix(SS[,"avC"],L,3); cSk[k,]=(SS[k,"sumC"]+(1:-1)*Si[1])/(SS[k,"nM"]+(1:-1));
			} else { cSk=matrix(grpS,L,3); cSk[k,]=(1-infl)*pSk[k,] + infl*(SS[k,"sumC"]+(1:-1)*Si[1])/SGC[k]; }
			groupPs=(A+diag(L))%*%apply(rbind(SGA,pSk,cSk), 2, function(x) payfun3(x[2:(L+1)],x[-(1:(L+1))],parms,pSG=x[1]))
			#Last step above adds group interaction effects into payoffs.
			#row entries in An are relative interaction strengths. Subtract diagonal to get mean Payoff difference
			An=t(t(A)%*%diag(1/(A%*%(1+0*A[,1]))[,1])); LambdsJ=((An-diag(L))%*%groupPs)[k,]; 
			dPg=c(dPg,parms["Alph"]*(LambdsJ[c(1,3)[1+(Mitigator)]]-LambdsJ[2]))
		}
	}
	freq=sum(SS[belonging,2])/sum(SS[belonging,1]); pOppose=c(freq,1-freq)[1+(Mitigator)] 
	#in next step, diff(Pays) = NM-M; want opposite when M is reconsidering
	Pr=parms["mu"] + (1-parms["mu"])*pOppose/(1+exp(parms["beta"]*((1-2*(Mitigator))*diff(Pays) + mean(dPg)))) 
	set.seed(seed); Swt=rbinom(1,1,prob=Pr); if(Swt==1) S[dSi,2]=S[dSi,1]*(!Mitigator); return(S);
}

#Run simulation
msim5=function(parms,Trun=1e2,seed0=15){
	ss0=function(S=seed0) set.seed(S); N=parms["N"]; L=parms["L"]; 
	ss0(); changes=c(1,cumsum(rbinom(Trun-1,N*L,prob=parms["xi"]))); Niter=tail(changes,1);
	ss0(); Seed=rnorm(Niter,sd=1e6); IDS=sample(1:(N*L),Niter,replace=TRUE); gID=floor((IDS-1)/N)+1; dS=cbind(gID,nID=IDS-(gID-1)*N);
	ss0(); X=mvrnorm(L*(L-1)/2, c(0,0), (1-parms["Acor"])*(diag(2)-1)+1); Ai=pnorm(X); 
	A=B=0*diag(L); A[upper.tri(A)]=Ai[,1]; B[upper.tri(B)]=Ai[,2]; A[lower.tri(A)]=t(B)[lower.tri(B)]; assign("A",A,.GlobalEnv);
	TS=group(parms,seed0); CS=sum(TS[,1]); out=sum(TS[,2])/CS;
	for(i in 2:Niter){ TS=m9iter(TS,parms,dSi=IDS[i],seed=Seed[i]); out=c(out,sum(TS[,2])/CS); }; 
	out=out[changes]; names(out)=changes; return(out);
}


#Example simulation with all processes: 
parms=Pmat[4,]; parms[c("Gcor","Glap","cvar","Alph","Glo","Tlocoff","infl")]=c(0.85,0.015,3,0.2,1,0.5,0.3); 
x=msim5(parms,Trun=150,seed0=20); plot(x,xlab="time step",ylab="global mitigation proportion")



#Running analyses above over ranges of parameters:
store=matrix(nrow=Trun,ncol=0);
for(i in 1:nrow(Starts)){
	print(i); pset=as.numeric(Starts[i,]); names(pset)=names(Starts); parms=Pmat[pset[1],];
	parms[c("T","cvar","Gcor","Glap","infl")]=c(pset[3:(ncol(Starts)-1)],0)[1:5]; 
	parms["c"]=pset[2]/(parms["rLt"]+(1-parms["rLt"])*(1-parms["rL"]))
	parms["Tlocoff"]=ifelse("Tlocoff"%in%names(pset),pset["Tlocoff"],parms["T"]);
	if("F"%in%names(pset)) parms["F"]=pset["F"]
	store=cbind(store, msim5(parms,Trun=Trun,seed0=1e4*pset["Rep"]))
}








