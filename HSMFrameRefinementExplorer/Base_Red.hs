module Base_Red where

-- Core part of the lambda phage, without reductions

----------------------------------------------
-- Initialization

emptypcond = [] :: [(Integer,Integer)]
initcond = [(1,emptypcond),(2,emptypcond)]

initstate = ([(0,0),(1,0),(2,0),(3,0)],initcond)

-- red, all states
sclist prot statevec = [(a,b)| a<-(bs prot), b<-(bs prot), (isHigh prot (lookval a statevec)) && not (isHigh prot (lookval b statevec))]

-- [(0,1),(2,0),(3,5)]

lookval a list = head [ v | (i,v) <- list, i == a]

isHigh 1 value = (value == 2) || (value == 3) || (value == 5)   -- cI
isHigh 2 value = (value == 4) || (value == 5) || (value == 3)   -- cro

statelist vlist = zip vlist (map statecond vlist)
statecond statevec = [(p,(sclist p statevec)) | p <- plistcond]  -- [1,2]

-- red
vectorlist = [[(0,a),(1,b),(2,c),(3,d)] | a<-[0,1], b<-[0,2,3,4,5],c<-[0,2,3,4,5],d<-[0,2,3,4,5] ]

----------------------------------------------

--red
bs 0 = [0]
bs 1 = [1,2,3]
bs 2 = [1,2,3]

plist = [0,1,2]    -- cII, cI, cro
plistcond = [1,2]

--red
pr s = if (bslookup 2 s == 0) && (bslookup 1 s == 0) then 1 else 0 
pl s = if (bslookup 4 s == 0) && (bslookup 5 s == 0) then 1 else 0

--bOR1[cI]*bOR3[] + bOR1[cro]*bOR3[] + bOR2[cI]*bOR3[] + bOR2[cro]*bOR3[] + bcII_1[cII]*bOR3[]
pme s = if (bslookup 3 s == 0) && ((bslookup 0 s == 1) || (bslookup 1 s) > 0 || (bslookup 2 s) > 0 ) then 1 else 0

-- [(0,a),(1,b),(2,c),(3,d)]
bslookup bs [] = -1
bslookup bs ((a,b):t) = if bs==a then b else bslookup bs t

pname a = case a of 
   0 -> "cII"
   1 -> "cI"
   2 -> "cro"
   3 -> "N"
   _ -> "??"