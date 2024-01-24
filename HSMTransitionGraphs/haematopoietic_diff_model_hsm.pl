

fname('h_diff_graph_hsm.txt').

% generate state id-s from given range

sid(2,X):-member(X,[0,1]).
sid(3,X):-member(X,[0,1,2]).

% generate transition from State to StateNext

trans(State,StateNext,1):-State=[Stat6,Il_4r,Il_4,Socs1,Gata3,Inf_y,Inf_yr,Stat1,T_bet],
						stat6(Stat6N,[Il_4r]),
						StateNext=[Stat6N,Il_4r,Il_4,Socs1,Gata3,Inf_y,Inf_yr,Stat1,T_bet].

trans(State,StateNext,2):-State=[Stat6,Il_4r,Il_4,Socs1,Gata3,Inf_y,Inf_yr,Stat1,T_bet],
						il_4r(Il_4rN,[Il_4,Socs1]),
						StateNext=[Stat6,Il_4rN,Il_4,Socs1,Gata3,Inf_y,Inf_yr,Stat1,T_bet].

trans(State,StateNext,3):-State=[Stat6,Il_4r,Il_4,Socs1,Gata3,Inf_y,Inf_yr,Stat1,T_bet],
						il_4(Il_4N,[Stat1,Gata3]),
						StateNext=[Stat6,Il_4r,Il_4N,Socs1,Gata3,Inf_y,Inf_yr,Stat1,T_bet].

trans(State,StateNext,4):-State=[Stat6,Il_4r,Il_4,Socs1,Gata3,Inf_y,Inf_yr,Stat1,T_bet],
						socs1(Socs1N,[Stat1,T_bet]),
						StateNext=[Stat6,Il_4r,Il_4,Socs1N,Gata3,Inf_y,Inf_yr,Stat1,T_bet].

trans(State,StateNext,5):-State=[Stat6,Il_4r,Il_4,Socs1,Gata3,Inf_y,Inf_yr,Stat1,T_bet],
						gata3(Gata3N,[Gata3,Stat6,T_bet]),
						StateNext=[Stat6,Il_4r,Il_4,Socs1,Gata3N,Inf_y,Inf_yr,Stat1,T_bet].

trans(State,StateNext,6):-State=[Stat6,Il_4r,Il_4,Socs1,Gata3,Inf_y,Inf_yr,Stat1,T_bet],
						inf_y(Inf_yN,[T_bet]),
						StateNext=[Stat6,Il_4r,Il_4,Socs1,Gata3,Inf_yN,Inf_yr,Stat1,T_bet].
						
trans(State,StateNext,7):-State=[Stat6,Il_4r,Il_4,Socs1,Gata3,Inf_y,Inf_yr,Stat1,T_bet],
						inf_yr(Inf_yrN,[Inf_y,Socs1]),
						StateNext=[Stat6,Il_4r,Il_4,Socs1,Gata3,Inf_y,Inf_yrN,Stat1,T_bet].

trans(State,StateNext,8):-State=[Stat6,Il_4r,Il_4,Socs1,Gata3,Inf_y,Inf_yr,Stat1,T_bet],
						stat1(Stat1N,[Inf_yr]),
						StateNext=[Stat6,Il_4r,Il_4,Socs1,Gata3,Inf_y,Inf_yr,Stat1N,T_bet].

trans(State,StateNext,9):-State=[Stat6,Il_4r,Il_4,Socs1,Gata3,Inf_y,Inf_yr,Stat1,T_bet],
						t_bet(T_betN,[T_bet,Stat1,Gata3]),
						StateNext=[Stat6,Il_4r,Il_4,Socs1,Gata3,Inf_y,Inf_yr,Stat1,T_betN].

% define arities of the involved genes

arities([2,2,2,2,2,3,3,3,3]).

% convert (somehow) state array to integer

s2numb(S,SNumb):-arities([A1,A2,A3,A4,A5,A6,A7,A8,A9]),S=[S1,S2,S3,S4,S5,S6,S7,S8,S9],
		SNumb is S1+A1*(S2+A2*(S3+A3*(S4+A4*(S5+A5*(S6+A6*(S7+A7*(S8+A8*(S9)))))))).
		
% write state - just output list of gene activities as a string
swrite([S1,S2,S3,S4,S5,S6,S7,S8,S9]):-write(S1),write(S2),write(S3),write(S4),write(S5),write(S6),write(S7),write(S8),write(S9).

% count the number of states for given arities vector
scount(Count,[A1,A2,A3,A4,A5,A6,A7,A8,A9]):-Count is A1*A2*A3*A4*A5*A6*A7*A8*A9.
	
% generate and write all state transitions

gen_graph:-arities([A1,A2,A3,A4,A5,A6,A7,A8,A9]),scount(Count,[A1,A2,A3,A4,A5,A6,A7,A8,A9]),write(Count),nl,
		sid(A1,S1),sid(A2,S2),sid(A3,S3),sid(A4,S4),sid(A5,S5),sid(A6,S6),sid(A7,S7),sid(A8,S8),sid(A9,S9),
		State=[S1,S2,S3,S4,S5,S6,S7,S8,S9],
		trans(State,StateNext,T),
		s2numb(State,StateNumb),s2numb(StateNext,StateNextNumb),
		write(StateNumb),tab(1),swrite(State),tab(1),write(StateNextNumb),tab(1),swrite(StateNext),
		%tab(1),write(T),
		nl,fail.
		

% just do it all :)
do_all:-fname(Fname),tell(Fname),gen_graph,told.		

	
		
% find state vector for a given state number (very inefficient, for manual checking only)

numb2s(SNumb):-arities([A1,A2,A3,A4,A5,A6,A7,A8,A9]),
		sid(A1,S1),sid(A2,S2),sid(A3,S3),sid(A4,S4),sid(A5,S5),sid(A6,S6),sid(A7,S7),sid(A8,S8),sid(A9,S9),
		State=[S1,S2,S3,S4,S5,S6,S7,S8,S9],
		s2numb(State,SNumb),writeln(State).
		






% stat6 = f[il_4r]
stat6(0,[0]).
stat6(1,[1]).

% il_4r = f[il_4,socs1]
il_4r(0,[0,0]).
il_4r(0,[0,1]).
il_4r(1,[1,0]).
il_4r(0,[1,1]).

% il_4 = f[stat1,gata3]
il_4(0,[0,0]).
il_4(1,[0,1]).
il_4(0,[1,0]).
il_4(0,[1,1]).
il_4(0,[2,0]).
il_4(0,[2,1]).

% socs1 = f[stat1,t_bet]
socs1(0,[0,0]).
socs1(1,[0,1]).
socs1(1,[0,2]).
socs1(1,[1,0]).
socs1(1,[1,1]).
socs1(1,[1,2]).
socs1(1,[2,0]).
socs1(1,[2,1]).
socs1(1,[2,2]).

% gata3 = f[gata3,stat6,t_bet]
gata3(0,[0,0,0]).
gata3(0,[0,0,1]).
gata3(0,[0,0,2]).
gata3(1,[0,1,0]).
gata3(0,[0,1,1]).
gata3(0,[0,1,2]).
gata3(1,[1,0,0]).
gata3(0,[1,0,1]).
gata3(0,[1,0,2]).
gata3(1,[1,1,0]).
gata3(0,[1,1,1]).
gata3(0,[1,1,2]).

% inf_y = f[t_bet]
inf_y(0,[0]).
inf_y(1,[1]).
inf_y(2,[2]).

% inf_yr = f[inf_y,socs1]
inf_yr(0,[0,0]).
inf_yr(0,[0,1]).
inf_yr(1,[1,0]).
inf_yr(1,[1,1]).
inf_yr(2,[2,0]).
inf_yr(1,[2,1]).

% stat1 = f[inf_yr]
stat1(0,[0]).
stat1(1,[1]).
stat1(2,[2]).

% t_bet = f[t_bet,stat1,gata3]
t_bet(0,[0,0,0]).
t_bet(0,[0,0,1]).
t_bet(1,[0,1,0]).
t_bet(0,[0,1,1]).
t_bet(2,[0,2,0]).
t_bet(0,[0,2,1]).
t_bet(1,[1,0,0]).
t_bet(0,[1,0,1]).
t_bet(1,[1,1,0]).
t_bet(0,[1,1,1]).
t_bet(2,[1,2,0]).
t_bet(0,[1,2,1]).
t_bet(2,[2,0,0]).
t_bet(0,[2,0,1]).
t_bet(2,[2,1,0]).
t_bet(0,[2,1,1]).
t_bet(2,[2,2,0]).
t_bet(0,[2,2,1]).
