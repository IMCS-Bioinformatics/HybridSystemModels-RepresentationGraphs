module Main where

import qualified Data.Map as Map
import qualified Data.Set as Set
import Data.Graph
import Data.Graph.SCC as SCC

-- Choice 1: choose the model
--import Base_RedEssential
--import Base_Red
import Base_BlueEssential
--import Base_Blue

-- Choice 2: choose the initialization
initlist = [initstate] -- zero-concentrations initial state
--initlist = (statelist vectorlist) -- any initial state

main = do 
   -- Choice 3: choose what to compute
   mapM_ print (xsort (do_scc cbb))   -- SCCs
   --mapM_ print (do_bisim bbb)         -- Bisimilarity classes
   --mapM_ print (xsort (do_states ss))   -- State space only

do_scc = condense_smapxc
do_bisim = condense_smapx
do_states = condense_statelist

condense_statevec [] = []
condense_statevec ((a,b):t) = (printone b) : (condense_statevec t)

clist ll = 
  [ (((pname a) ++ (chdescr b)),c,((condense_statevec d)++(changevec d),(compact_cond_block e))) | (a,b,c,(d,e)) <- ll ]

printone x = case x of 
   0 -> '0'
   1 -> '1'
   2 -> '2'
   3 -> '3'
   4 -> '4'
   5 -> '5'
   6 -> '6'   -- not needed right now
   _ -> '?'

cstate ((a,c),b) = ((condense_statevec a)++(changevec a),(compact_cond_block c),(clist b))
condense_statelist ss = [ (cstate x) | x <- ss ]

cstate_short ((a,c),b) = ((condense_statevec a),(compact_cond_block c))
cstate_short_plain (a,c) = ((condense_statevec a)++(changevec a),(compact_cond_block c))

-- compact_cond_block at the end of the file (build the compact linearized form)

--condense_smap_item (a,blist) = ((cstate a), (map cstate_short blist))
condense_smap_item (a,blist) = (cstate a)

condense_smap = map condense_smap_item

condense_smapx_item (a,slist,linklist) = 
   ((cstate_short_plain a), (map cstate_short_plain slist), (clist linklist))
condense_smapx = map condense_smapx_item

condense_smapxc_item (c, a,slist,linklist) = 
   (c,(cstate_short_plain a), (map cstate_short_plain slist), (clist linklist))
condense_smapxc = map condense_smapxc_item
----------------------------------------------
-- Build bisimulation

-- Initial markup of candidates
bisimap_start = 
   let ssi = xsort ss in
     [(a,b,1) | a <- ssi, b <- ssi, (statesmatch a b) ]

statesmatch a b = (basestate (fst (statepart a))) == (basestate (fst (statepart b)))

bisimap_next bmap = 
   [(a,b,(if x==0 then 0 else if (match a b bmap) then 1 else 0))| (a,b,x) <- bmap ] 

match (ssa,alist) (ssb,blist) bmap = (dmatch ssa alist blist bmap) && (dmatch ssb blist alist bmap)

dmatch ssa [] blist bmap = True
dmatch ssa (move : rest) blist bmap = 
   if (not ((basestate (fst ssa)) == (basestate (fst (nodepart move))))) && (not (ematch move blist bmap)) 
        then False else dmatch ssa rest blist bmap
   -- remove self-moves from matching requirement (does not help, in fact)

ematch (p,ch,bs,state) blist bmap = 
   let bmatches = [ state1 | (p1,ch1,bs1,state1) <- blist, p1==p && ch1==ch && bs1==bs1 ] in
     if bmatches == [] then False
       else let realmatches = [ (a,b,x) | (a,b,x) <- bmap, state1 <- bmatches, 
                                   statepart a == state && statepart b == state1 && x == 1 ] in
           if realmatches == [] then False else True

matchcount [] = 0
matchcount ((_,_,x) : tail) =  x + (matchcount tail)

bisilist = bisimap_start : (map bisimap_next bisilist)   -- Infinite iteration list

bisimap_find (a:b:rest) = if matchcount(a)==matchcount(b) then a else bisimap_find (b:rest)

-- bisimap = xsort [(a,b) | (a,b,c) <- (bisimap_find bisilist), (c==1 && a <= b) ]
bisimap = xsort [(a,b) | (a,b,c) <- (bisimap_find bisilist), c==1]

-----------------------------------------------
-- process bisimap to produce equivalence classes
sortAndGroup assocs = Map.fromListWith (++) [(k, [v]) | (k, v) <- assocs]

bb = [ (a,bstates) | (a,bstates) <- (Map.toAscList (sortAndGroup bisimap)), (bgroup_ok a bstates)]

bgroup_ok a [] = True
bgroup_ok a (b:t) = if a > b then False else (bgroup_ok a t)

-----------------------------------------------
-- build links between equivalence classes

-- pointers to the representatives of the equivalence class
revpair (a,b) = [((fst x),(fst a)) | x <- b]
revlist [] = []
revlist (h:t) = (revpair h) ++ (revlist t)

collect_links [] = []
collect_links (h:t) = (snd h) ++ (collect_links t)

update_links linkcoll rlist = [ (p,ch,bs,b) | (p,ch,bs,ss) <- linkcoll, (a,b) <- rlist, ss==a]

grouplinks blist rlist = (dedup (xsort (update_links (collect_links blist) rlist)))

bbb = let rlist = revlist bb in
  [ (state, (map fst bstates), (grouplinks bstates rlist))| ((state,transitions),bstates) <- bb ]
 -- Take (a,blist) from bb
 -- Return a, the other sites and updated links. rlist is the dictionary for list update.

----------------------------------------------
-- Strongly connected components
ge = [let num = getnum state; numlinks = getnumlist links in (num,num,numlinks) | (state, otherstates, links) <- bbb ]
teststate = let (a,b,links) = head bbb in links

getnum state = Map.findWithDefault 0 state ss_to_num
getnumlist links = map getnum (map nodepart links)

(graph,nodeFromVertex,vertexFromKey) = graphFromEdges ge

getkey state = let nn = (getnum state) in case (vertexFromKey nn) of {Just x -> x; _ -> 0 }

(comp_vlist,gcomponent) = SCC.scc graph

scomponent state = gcomponent (getkey state)
link_to_component (a,b,c,d) = (a,b,c,(scomponent d))

cbb = [((scomponent state), state, otherstates, links) | (state, otherstates, links) <- bbb ]

cbx = [ (comp,(state,otherstates,(map link_to_component links))) | (comp,state,otherstates,links) <- cbb ]

component_map = sortAndGroup cbx

component_list = Map.toAscList component_map

mergelists elist = dedup (foldr (++) [] [linklist |(state,otherstates,linklist) <- elist ])

xclist = [(comp, (mergelists elist), elist) | (comp,elist) <- component_list]

xcsimplelist = [(a,b) | (a,b,c) <- xclist]

-----------------------------------------------
-- generate state space with conditions and transitions

natnums = [1..]
ss_to_num = Map.fromList (zip (fst sx) natnums)
num_to_ss = Map.fromList (zip natnums (fst sx))

sx = let (sl,mm) = stateset Set.empty Map.empty (Set.fromList initlist) in ((Set.toList sl),(Map.toList mm))
ss = snd sx
--ss = let (sl,mm) = stateset Set.empty Map.empty (Set.empty) in (Map.toList mm)

stateset hi li addons = 
  if (Set.null addons) then (hi,li)
  else let (h,tail) = ((Set.elemAt 0 addons),(Set.deleteAt 0 addons)) in
    let nn = (allmoves h) in
      let nx = dedup [ x | x <- (map nodepart nn), 
           (Set.notMember x hi) && (Set.notMember x addons)] in
        stateset (Set.insert h hi) (Map.insert h nn li) (Set.union (Set.fromList nx) tail)

nodepart (_,_,_,x) = x    -- Gets node description out of (p,ch,bs,node)
statepart (x,_) = x       -- Gets state description out of (state,transitions)

enrich ss nn = (ss,nn)

-- All moves from a state, by all proteins
allmoves ss = foldr (++) [] [consistent_moves ss p | p <- plist]

update_state [] p ch bb = []
update_state ((bs,mode):tail) p ch bb = 
  if bb==bs 
    then (bs,(nextmode p ch mode)) : tail 
    else (bs,mode) : (update_state tail p ch bb)

-- Move from a state by a protein and binding site
simple_move (s,c) p bs = let ch = (change s p) in (p, ch, bs, ((update_state s p ch bs),c))

-- Make list of all moves, exclude self-moves
simple_moves ss p = let movelist = [ simple_move ss p bs | bs <- (bs p)] 
                        in [ x | x <- movelist, not (nodepart x == ss) ]

-- Moves with plain merge of conditions
moves ss p = let mlist = simple_moves ss p in 
  let bslist = [ bb | (_,_,bb,_) <- mlist ] in
    [(p,ch,b,(s,(update_condlist sclist p ch b bslist))) | (p,ch,b,(s,sclist)) <- mlist]
--[(p,ch,b,(s,[(cmerge c (condlist p ch b bslist))| (pp,c) <- sclist ])) | (p,ch,b,(s,sclist)) <- mlist]

update_condlist state_condlist p ch b bslist = 
  [ (pp, if pp==p then (cmerge x (condlist ch b bslist)) else x ) | (pp,x) <- state_condlist ]

-- conditions comparing this move with other possible moves by this protein at other binding sites
condlist ch bs bslist = [ buildcond ch bs b  | b <- bslist, not (b == bs) ]

-- build a condition comparing this move with the other possible move
buildcond ch bs b = if ch==1 then (bs,b) else (b,bs)

cmerge a b = let x = (xsort (dedup (a ++ b))) in xsort (ifind (ii x))
--cmerge a b = let x = (xsort (dedup (a ++ b))) in checkloop(ifind (ii x))

ii x = x : (map (inext x) (ii x))
inext x y = dedup (y ++ (compose y x))
ifind (a:b:list) = if (length a)==(length b) then a else ifind (b:list)

compose x y = [(a,d) | (a,b) <- x, (c,d) <- y, b==c]

hasloop [] = False
hasloop ((a,b):t) = if a==b then True else hasloop t

-- checkloop x = if (hasloop x) then [(-1,-1)] else x		-- Error value introduced

consistent_moves ss p = [(p,ch,b,(s,c)) | (p,ch,b,(s,c)) <- (moves ss p), (consistent c)]

consistent [] = True
consistent ((a,b):t) = if (hasloop b) then False else (consistent t)

-- iconsistent clist = hasloop c

-- iconsistent c = rconsistent (c ++ (compose c (c ++ (compose c c)) ))
-- compose x y = [(a,d) | (a,b) <- x, (c,d) <- y, b==c]

-- rconsistent [] = True
-- rconsistent ((a,b):t) = if a==b then False else rconsistent t

change s p = case p of   -- growth/degradation of a protein in a BS
    1 -> pme s           -- cI				
    3 -> pl s            -- N
    _ -> pr s            -- cro, cII

chdescr 0 = "-"
chdescr 1 = "+"

changevec state = foldr (++) [] (map chdescr (map (change state) plist))

-- moves of a binding site by a protein value change
-- for multi-protein BS: 0 - empty, 2 - cI, 4 - cro, 3 - cI + hidden cro, 5 - cro + hidden cI
nextmode p ch mode = case (p,ch) of 
    (1,0) -> case mode of 2 -> 0
                          3 -> 4        -- if cro was behind cI
                          _ -> mode
    (1,1) -> case mode of 0 -> 2
                          4 -> 5        -- cI+, hidden behind cro
                          _ -> mode
    (2,0) -> case mode of 4 -> 0
                          5 -> 2        -- if cI was behind cro
                          _ -> mode
    (2,1) -> case mode of 0 -> 4
                          2 -> 3        -- cro+, hidden behind cI
                          _ -> mode
    _ -> ch

basemode n = case n of 3 -> 2 
                       5 -> 4
                       _ -> n

basestate [] = []
basestate ((a,b):t) = (a, (basemode b)) : (basestate t)

-- auxiliaries
inlist x [] = False
inlist x (h:t) = if x==h then True else (inlist x t)

dedup [] = []
dedup (h:t) = if (inlist h t) then (dedup t) else h:(dedup t)

psort [] = []
psort (h:t) = pinsert h (psort t)

pinsert (a,b) [] = [(a,b)]
pinsert (a,b) ((c,d):tt) = if a < c || (a == c && b < d) then (a,b) : ((c,d):tt) else (c,d) : (pinsert (a,b) tt)

xsort [] = []
xsort (h:t) = xinsert h (xsort t)

xinsert a [] = [a]
xinsert a (h:t) = if a < h then a:(h:t) else h:(xinsert a t)

----------------------------------------------------
--compact condition presentation 
compact_cond_block [] = []
compact_cond_block ((a,clist):t) = 
   (stringify_cond_list ((pname a),(compact_cond_list clist))) : (compact_cond_block t)

compact_cond_list ll = reverse (c_cond_list (childsort ll ll))

stringify [] = []
stringify (h:t) = (printone h) : (stringify t) 

stringify_list [] = []
stringify_list (h:[]) = (stringify h)
stringify_list (h:t) = (stringify h) ++ "." ++ (stringify_list t)

stringify_cond_list (a,b) = a ++ ":" ++ stringify_list b

c_cond_list [] = [] 
c_cond_list ((a,b):tail) = buildchainlist [] [a,b] b tail

-- sorting (a,b) descending by the number of descendants (edges (a,x)) in ll for a, then for b
childsort [] ll = []
childsort (h:t) ll = ch_insert ll h (childsort t ll)

ch_insert ll (a,b) [] = (a,b):[]
ch_insert ll (a,b) ((c,d):tt) = 
   let {ca = (childcount a ll); cb = (childcount b ll); cc = (childcount c ll); cd = (childcount d ll)} in
     let bb = if ca > cc || (ca == cc && cb > cd) then True 
                 else if ca==cc && cb==cd && (a < c || (a==c && b < d)) then True else False in
      if bb then (a,b):(c,d):tt else (c,d):(ch_insert ll (a,b) tt)

childcount c [] = 0
childcount c ((a,b):t) = let n = (childcount c t) in if c==a then n+1 else n

buildchainlist chainlist chain _ [] = chain:chainlist

buildchainlist chainlist chain next list = 
   let (a,b) = xlookup next list in  -- find the next tuple to be added 
      if next == a then              -- the current chain can be continued 
        -- Add b to the chain, remove all (x,b) for x in the chain, from the list
        (buildchainlist chainlist (chain++[b]) b (clearlist chain b list))
      else (buildchainlist (chain:chainlist) [a,b] b (tail list))  
           -- move the chain to chainlist, start a new chain

xlookup x pairlist = 
   let ll = [(a,b) | (a,b) <- pairlist, a == x] in
     if ll == [] then (head pairlist) else (head ll)

clearlist chain next list = 
   [(a,b) | (a,b) <- list, (b /= next || (notin a (next:chain)) )]

notin a [] = True
notin a (h:t) = if a==h then False else (notin a t)
----------------------------------------------