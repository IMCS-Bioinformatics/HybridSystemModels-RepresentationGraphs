fname('h_diff_graph_simple.txt').
fname_a('h_diff_graph_hsm.txt').
fname_a2('h_diff_graph_complete.txt').
fname_rgb('h_diff_hsm_base_sg.txt').
fname_rgv('h_diff_hsm_variant_sg.txt').

% heamatopoietic cell differentiation model

% generate state id-s from given range

sid(2,X):-member(X,[0,1]).
sid(3,X):-member(X,[0,1,2]).

% generate transition from State to StateNext

trans(State,StateNext):-State=[Gata_2,Gata_1,C_ebpa,Pu_1,Egrnab,Eklf,Fli_1,Scl,Cjun,Gfi_1,Fog_1],
						gata_2(Gata_2N,[Gata_2,Gata_1,Fog_1,Pu_1]),
						gata_1(Gata_1N,[Gata_1,Gata_2,Fli_1,Pu_1]),
						c_ebpa(C_ebpaN,[C_ebpa,Gata_1,Fog_1,Scl]),
						pu_1(Pu_1N,[C_ebpa,Pu_1,Gata_1,Gata_2]),
						egrnab(EgrnabN,[Pu_1,Cjun,Gfi_1]),
						eklf(EklfN,[Gata_1,Fli_1]),
						fli_1(Fli_1N,[Gata_1,Eklf]),
						scl(SclN,[Gata_1,Pu_1]),
						cjun(CjunN,[Pu_1,Gfi_1]),
						gfi_1(Gfi_1N,[C_ebpa,Egrnab]),
						fog_1(Fog_1N,[Gata_1]),
						StateNext=[Gata_2N,Gata_1N,C_ebpaN,Pu_1N,EgrnabN,EklfN,Fli_1N,SclN,CjunN,Gfi_1N,Fog_1N].

% generate transition from State to StateNext applying ANY SUBSET of rules

trans_set(State,StateNext):-State=[Gata_2,Gata_1,C_ebpa,Pu_1,Egrnab,Eklf,Fli_1,Scl,Cjun,Gfi_1,Fog_1],
						(gata_2(Gata_2N,[Gata_2,Gata_1,Fog_1,Pu_1]);Gata_2N=Gata_2),
						(gata_1(Gata_1N,[Gata_1,Gata_2,Fli_1,Pu_1]);Gata_1N=Gata_1),
						(c_ebpa(C_ebpaN,[C_ebpa,Gata_1,Fog_1,Scl]);C_ebpaN=C_ebpa),
						(pu_1(Pu_1N,[C_ebpa,Pu_1,Gata_1,Gata_2]);Pu_1N=Pu_1),
						(egrnab(EgrnabN,[Pu_1,Cjun,Gfi_1]);EgrnabN=Egrnab),
						(eklf(EklfN,[Gata_1,Fli_1]);EklfN=Eklf),
						(fli_1(Fli_1N,[Gata_1,Eklf]);Fli_1N=Fli_1),
						(scl(SclN,[Gata_1,Pu_1]);SclN=Scl),
						(cjun(CjunN,[Pu_1,Gfi_1]);CjunN=Cjun),
						(gfi_1(Gfi_1N,[C_ebpa,Egrnab]);Gfi_1N=Gfi_1),
						(fog_1(Fog_1N,[Gata_1]);Fog_1N=Fog_1),
						StateNext=[Gata_2N,Gata_1N,C_ebpaN,Pu_1N,EgrnabN,EklfN,Fli_1N,SclN,CjunN,Gfi_1N,Fog_1N].
						
% just in case rules for asynchronous generation					

trans(State,StateNext,1,Gata_2N):-State=[Gata_2,Gata_1,C_ebpa,Pu_1,Egrnab,Eklf,Fli_1,Scl,Cjun,Gfi_1,Fog_1],
						gata_2(Gata_2N,[Gata_2,Gata_1,Fog_1,Pu_1]),
						StateNext=[Gata_2N,Gata_1,C_ebpa,Pu_1,Egrnab,Eklf,Fli_1,Scl,Cjun,Gfi_1,Fog_1].

trans(State,StateNext,2,Gata_1N):-State=[Gata_2,Gata_1,C_ebpa,Pu_1,Egrnab,Eklf,Fli_1,Scl,Cjun,Gfi_1,Fog_1],
						gata_1(Gata_1N,[Gata_1,Gata_2,Fli_1,Pu_1]),
						StateNext=[Gata_2,Gata_1N,C_ebpa,Pu_1,Egrnab,Eklf,Fli_1,Scl,Cjun,Gfi_1,Fog_1].
						
trans(State,StateNext,4,C_ebpaN):-State=[Gata_2,Gata_1,C_ebpa,Pu_1,Egrnab,Eklf,Fli_1,Scl,Cjun,Gfi_1,Fog_1],
						c_ebpa(C_ebpaN,[C_ebpa,Gata_1,Fog_1,Scl]),
						StateNext=[Gata_2,Gata_1,C_ebpaN,Pu_1,Egrnab,Eklf,Fli_1,Scl,Cjun,Gfi_1,Fog_1].

trans(State,StateNext,8,Pu_1N):-State=[Gata_2,Gata_1,C_ebpa,Pu_1,Egrnab,Eklf,Fli_1,Scl,Cjun,Gfi_1,Fog_1],
						pu_1(Pu_1N,[C_ebpa,Pu_1,Gata_1,Gata_2]),
						StateNext=[Gata_2,Gata_1,C_ebpa,Pu_1N,Egrnab,Eklf,Fli_1,Scl,Cjun,Gfi_1,Fog_1].
						
trans(State,StateNext,16,EgrnabN):-State=[Gata_2,Gata_1,C_ebpa,Pu_1,Egrnab,Eklf,Fli_1,Scl,Cjun,Gfi_1,Fog_1],
						egrnab(EgrnabN,[Pu_1,Cjun,Gfi_1]),
						StateNext=[Gata_2,Gata_1,C_ebpa,Pu_1,EgrnabN,Eklf,Fli_1,Scl,Cjun,Gfi_1,Fog_1].

trans(State,StateNext,32,EklfN):-State=[Gata_2,Gata_1,C_ebpa,Pu_1,Egrnab,Eklf,Fli_1,Scl,Cjun,Gfi_1,Fog_1],
						eklf(EklfN,[Gata_1,Fli_1]),
						StateNext=[Gata_2,Gata_1,C_ebpa,Pu_1,Egrnab,EklfN,Fli_1,Scl,Cjun,Gfi_1,Fog_1].
						
trans(State,StateNext,64,Fli_1N):-State=[Gata_2,Gata_1,C_ebpa,Pu_1,Egrnab,Eklf,Fli_1,Scl,Cjun,Gfi_1,Fog_1],
						fli_1(Fli_1N,[Gata_1,Eklf]),
						StateNext=[Gata_2,Gata_1,C_ebpa,Pu_1,Egrnab,Eklf,Fli_1N,Scl,Cjun,Gfi_1,Fog_1].						

trans(State,StateNext,128,SclN):-State=[Gata_2,Gata_1,C_ebpa,Pu_1,Egrnab,Eklf,Fli_1,Scl,Cjun,Gfi_1,Fog_1],
						scl(SclN,[Gata_1,Pu_1]),
						StateNext=[Gata_2,Gata_1,C_ebpa,Pu_1,Egrnab,Eklf,Fli_1,SclN,Cjun,Gfi_1,Fog_1].
						
trans(State,StateNext,256,CjunN):-State=[Gata_2,Gata_1,C_ebpa,Pu_1,Egrnab,Eklf,Fli_1,Scl,Cjun,Gfi_1,Fog_1],
						cjun(CjunN,[Pu_1,Gfi_1]),
						StateNext=[Gata_2,Gata_1,C_ebpa,Pu_1,Egrnab,Eklf,Fli_1,Scl,CjunN,Gfi_1,Fog_1].

trans(State,StateNext,512,Gfi_1N):-State=[Gata_2,Gata_1,C_ebpa,Pu_1,Egrnab,Eklf,Fli_1,Scl,Cjun,Gfi_1,Fog_1],
						gfi_1(Gfi_1N,[C_ebpa,Egrnab]),
						StateNext=[Gata_2,Gata_1,C_ebpa,Pu_1,Egrnab,Eklf,Fli_1,Scl,Cjun,Gfi_1N,Fog_1].
						
trans(State,StateNext,1024,Fog_1N):-State=[Gata_2,Gata_1,C_ebpa,Pu_1,Egrnab,Eklf,Fli_1,Scl,Cjun,Gfi_1,Fog_1],
						fog_1(Fog_1N,[Gata_1]),
						StateNext=[Gata_2,Gata_1,C_ebpa,Pu_1,Egrnab,Eklf,Fli_1,Scl,Cjun,Gfi_1,Fog_1N].	
			
% specify which rules will be used
rule_list([1,2,4,8,16,32,64,128,256,512,1024]).
			
% define arities of the involved genes

arities([2,2,2,2,2,2,2,2,2,2,2]).

% convert (somehow) state array to integer

s2numb(S,SNumb):-arities([A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11]),S=[S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11],
		SNumb is S11+A11*(S10+A10*(S9+A9*(S8+A8*(S7+A7*(S6+A6*(S5+A5*(S4+A4*(S3+A3*(S2+A2*(S1)))))))))).
		
% write state - just output list of gene activities as a string
swrite([S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11]):-write(S1),write(S2),write(S3),write(S4),write(S5),write(S6),write(S7),write(S8),write(S9),write(S10),write(S11).

% count the number of states for given arities vector
scount(Count,[A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11]):-Count is A1*A2*A3*A4*A5*A6*A7*A8*A9*A10*A11.
	
% generate and write all state transitions

gen_graph:-arities([A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11]),scount(Count,[A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11]),write(Count),nl,
		sid(A1,S1),sid(A2,S2),sid(A3,S3),sid(A4,S4),sid(A5,S5),sid(A6,S6),sid(A7,S7),sid(A8,S8),sid(A9,S9),sid(A10,S10),sid(A11,S11),
		State=[S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11],
		trans(State,StateNext), % normal complete version
		s2numb(State,StateNumb),s2numb(StateNext,StateNextNumb),
		write(StateNumb),tab(1),swrite(State),tab(1),write(StateNextNumb),tab(1),swrite(StateNext),tab(1),write('0'),nl,fail.
gen_graph.
		
% asyncronous version assuming only a single rulle trigered transitions from states
		
gen_graph_a:-arities([A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11]),scount(Count,[A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11]),write(Count),nl,
		sid(A1,S1),sid(A2,S2),sid(A3,S3),sid(A4,S4),sid(A5,S5),sid(A6,S6),sid(A7,S7),sid(A8,S8),sid(A9,S9),sid(A10,S10),sid(A11,S11),
		State=[S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11],
		rule_list(RL),
		trans(State,StateNext,T,V), % asyncrhronous version with independent transitions
		member(T,RL),
		%T \= k,
		switch(T,ST,V),
		s2numb(State,StateNumb),s2numb(StateNext,StateNextNumb),
		write(StateNumb),tab(1),swrite(State),tab(1),write(StateNextNumb),tab(1),swrite(StateNext),tab(1),write(ST),nl,fail.
gen_graph_a.

% switch T value sign to negative for transitions 1->0

switch(T,T,1).
switch(T,ST,0):-ST is -T.

% get gene nummber from rule (admittedly not the msot efficient way to do this)

gene_numb(1,0).
gene_numb(2,1).
gene_numb(4,2).
gene_numb(8,3).
gene_numb(16,4).
gene_numb(32,5).
gene_numb(64,6).
gene_numb(128,7).
gene_numb(256,8).
gene_numb(512,9).
gene_numb(1024,10).

trans_sign(1,'+').
trans_sign(0,'-').

% single rule asyncronous version with return of values instead of printing them

gen_graph_a_ret([State,StateNumb,StateNext,StateNextNumb,GN,TS]):-arities([A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11]),
		sid(A1,S1),sid(A2,S2),sid(A3,S3),sid(A4,S4),sid(A5,S5),sid(A6,S6),sid(A7,S7),sid(A8,S8),sid(A9,S9),sid(A10,S10),sid(A11,S11),
		State=[S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11],
		rule_list(RL),
		trans(State,StateNext,T,V), % asyncrhronous version with independent transitions
		member(T,RL),
		%T \= k,
		%switch(T,ST,V),
		% get a gene # and transition sign +/- instead
		gene_numb(T,GN),
		trans_sign(V,TS),
		s2numb(State,StateNumb),s2numb(StateNext,StateNextNumb).

% asyncronous version assuming any subset of rules

gen_graph_a2:-arities([A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11]),scount(Count,[A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11]),write(Count),nl,
		sid(A1,S1),sid(A2,S2),sid(A3,S3),sid(A4,S4),sid(A5,S5),sid(A6,S6),sid(A7,S7),sid(A8,S8),sid(A9,S9),sid(A10,S10),sid(A11,S11),
		State=[S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11],
		trans_set(State,StateNext), % apply all possible subsets of transitions (including empty set)
		s2numb(State,StateNumb),s2numb(StateNext,StateNextNumb),
		write(StateNumb),tab(1),swrite(State),tab(1),write(StateNextNumb),tab(1),swrite(StateNext),tab(1),write('0'),nl,fail.
gen_graph_a2.		
		
% just do it all :)
do_all:-fname(Fname),tell(Fname),gen_graph,told.	

% asynchronous version	
do_all_a:-fname_a(Fname),tell(Fname),gen_graph_a,told.

% asynchronous version using all subsets	
do_all_a2:-fname_a2(Fname),tell(Fname),gen_graph_a2,told.	

% write the output in syntax compatible with representation graph generation input format
do_rg:-fname_rgb(Fname),tell(Fname),bagof(X,gen_graph_a_ret(X),L),length(L,Len),write_rg_header(Len),write_rg_list(L),told.   	

% write file header in representation representation graph generation input format
write_rg_header(Len):-
	writeln('### Myeloid differentiation model state graph - base version ###'),
	%writeln('### Myeloid differentiation model state graph - variant without direct EgrNab inhibition ###'),
	writeln('#'),
	writeln('# number of genes | int'),
	writeln(11), % currently constant
	writeln('# number of binding sites | int'),
	writeln(11), % currently constant
	writeln('# number of states | int'), 
	arities([A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11]),scount(Count,[A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11]),
	writeln(Count),
	writeln('# 0 state id | int'), % define 0 as initial state
	writeln(0),
	writeln('# gene names | str[]'),
	writeln('GATA-2 GATA-1 C/EBPa PU.1 EgrNab EKLF Fli-1 SCL cJun Gfi-1 FOG-1'), % currently constant
	writeln('# binding site names | str[]'),
	writeln('bGATA-2 bGATA-1 bC/EBPa bPU.1 bEgrNab bEKLF bFli-1 bSCL bcJun bGfi-1 bFOG-1'),nl, % currently constant
	writeln('# list of states and transitions'), 
	writeln('# state attributes: id gene_states bs_states type edge_n | int str str int int'),
	writeln('# transition attributes: dest_id gene_id [+|-] | int int def'),nl.

% write transition list in representation representation graph generation input format
% currently hardcoded assuming 11 sequential transition batches for each of the states

write_rg_list([]).
write_rg_list([[S,Sn,NS00,NSn00,GN00,TS00],[S,Sn,NS01,NSn01,GN01,TS01],[S,Sn,NS02,NSn02,GN02,TS02],[S,Sn,NS03,NSn03,GN03,TS03],
	[S,Sn,NS04,NSn04,GN04,TS04],[S,Sn,NS05,NSn05,GN05,TS05],[S,Sn,NS06,NSn06,GN06,TS06],[S,Sn,NS07,NSn07,GN07,TS07],
	[S,Sn,NS08,NSn08,GN08,TS08],[S,Sn,NS09,NSn09,GN09,TS09],[S,Sn,NS10,NSn10,GN10,TS10]|Tail]):-
	write(Sn),tab(1),write_list(S),tab(1),write_list(S),tab(1),state_type(Sn,St),tab(1),write(St),tab(1),write(11),nl, % assumes here 11 transitions will follow
    write(NSn00),tab(1),write(GN00),tab(1),write(TS00),tab(1),	
    write(NSn01),tab(1),write(GN01),tab(1),write(TS01),tab(1),	
    write(NSn02),tab(1),write(GN02),tab(1),write(TS02),tab(1),	
    write(NSn03),tab(1),write(GN03),tab(1),write(TS03),tab(1),	
    write(NSn04),tab(1),write(GN04),tab(1),write(TS04),tab(1),	
    write(NSn05),tab(1),write(GN05),tab(1),write(TS05),tab(1),	
    write(NSn06),tab(1),write(GN06),tab(1),write(TS06),tab(1),	
    write(NSn07),tab(1),write(GN07),tab(1),write(TS07),tab(1),	
    write(NSn08),tab(1),write(GN08),tab(1),write(TS08),tab(1),	
    write(NSn09),tab(1),write(GN09),tab(1),write(TS09),tab(1),	
    write(NSn10),tab(1),write(GN10),tab(1),write(TS10),nl,nl,
	write_rg_list(Tail).		
	
% write list elements in a sequential string
write_list([]).
write_list([H|T]):-write(H),write_list(T).

% state type - retrurns 1 for known steady states and 0 otherwise
state_type(386,1):-!. % Gran
state_type(452,1):-!. % Mono
state_type(537,1):-!. % Mega
state_type(553,1):-!. % Ery
state_type(_,0). % other
		
% find state vector for a given state number (very inefficient, for manual checking only)

numb2s(SNumb):-arities([A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11]),
		sid(A1,S1),sid(A2,S2),sid(A3,S3),sid(A4,S4),sid(A5,S5),sid(A6,S6),sid(A7,S7),sid(A8,S8),sid(A9,S9),sid(A10,S10),sid(A11,S11),
		State=[S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11],
		s2numb(State,SNumb),writeln(State).
		

% try to define some logical operators: && (and), ++ (or), \\ (not) and eq for their computation

eq(0,0):-!.
eq(1,1):-!.
eq(0,\\(1)):-!.
eq(1,\\(0)):-!.
eq(R,\\(X)):-eq(X1,X),eq(R,\\(X1)).
eq(0,&&(0,0)):-!.
eq(0,&&(0,1)):-!.
eq(0,&&(1,0)):-!.
eq(1,&&(1,1)):-!.
eq(R,&&(X,Y)):-eq(X1,X),eq(Y1,Y),eq(R,&&(X1,Y1)).
eq(0,++(0,0)):-!.
eq(1,++(0,1)):-!.
eq(1,++(1,0)):-!.
eq(1,++(1,1)):-!.
eq(R,++(X,Y)):-eq(X1,X),eq(Y1,Y),eq(R,++(X1,Y1)).
% define operator syntax for these
:-op(700,xfx,eq).
:-op(600,fy,\\).
:-op(500,yfx,++).
:-op(400,yfx,&&).



%Gata_2 = f[Gata_2,Gata_1,Fog_1,Pu_1]
gata_2(R,[Gata_2,Gata_1,Fog_1,Pu_1]):-R eq Gata_2 && \\(Gata_1 && Fog_1) && \\(Pu_1).

%Gata_1 = f[Gata_1,Gata_2,Fli_1,Pu_1]
gata_1(R,[Gata_1,Gata_2,Fli_1,Pu_1]):-R eq (Gata_1 ++ Gata_2 ++ Fli_1) && \\(Pu_1).

%C_ebpa = f[C_ebpa,Gata_1,Fog_1,Scl]
c_ebpa(R,[C_ebpa,Gata_1,Fog_1,Scl]):-R eq C_ebpa && \\(Gata_1 && Fog_1 && Scl).

%Pu_1 = f[C_ebpa,Pu_1,Gata_1,Gata_2]
pu_1(R,[C_ebpa,Pu_1,Gata_1,Gata_2]):-R eq (C_ebpa ++ Pu_1) && \\(Gata_1 ++ Gata_2).

%Egrnab = f[Pu_1,Cjun,Gfi_1]
% base assumption
egrnab(R,[Pu_1,Cjun,Gfi_1]):-R eq (Pu_1 && Cjun) && \\(Gfi_1).
% variant without direct inhibition
%egrnab(R,[Pu_1,Cjun,Gfi_1]):-R eq (Pu_1 && Cjun).

%Eklf = f[Gata_1,Fli_1]
eklf(R,[Gata_1,Fli_1]):-R eq Gata_1 && \\(Fli_1).

%Fli_1 = f[Gata_1,Eklf]
fli_1(R,[Gata_1,Eklf]):-R eq Gata_1 && \\(Eklf).

%Scl = f[Gata_1,Pu_1]
scl(R,[Gata_1,Pu_1]):-R eq Gata_1 && \\(Pu_1).

%Cjun = f[Pu_1,Gfi_1]
cjun(R,[Pu_1,Gfi_1]):-R eq Pu_1 && \\(Gfi_1).

%Gfi_1 = f[C_ebpa,Egrnab]
gfi_1(R,[C_ebpa,Egrnab]):-R eq  C_ebpa && \\(Egrnab).

%Fog_1 = f[Gata_1]
fog_1(R,[Gata_1]):-R eq Gata_1.

