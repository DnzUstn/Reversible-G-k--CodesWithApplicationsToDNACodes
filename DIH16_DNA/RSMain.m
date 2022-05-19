// *************************************************************************************************
// ***** Copyright(c) 2022, Tarsus University, Department of Mathematical and Computer Science *****
// *************************************************************************************************
//  Project Code  : PRJ_RS-001
//  Project Title : DIH-16 DNA Codes 
//  Publisher     : Assoc.Prof.Dr. Serap Şahinkaya, Assoc.Prof.Dr. Deniz Ustun and Dr. Adrian Korban 
//  Developer     : Assoc.Prof.Dr. Serap Şahinkaya, Assoc.Prof.Dr. Deniz Ustun and Dr. Adrian Korban
//  Contact Info  : Assoc.Prof.Dr. Serap Şahinkaya, serap@tarsus.edu.tr
//                  Assoc.Prof.Dr. Deniz Ustun, denizustun@tarsus.edu.tr
//                             Dr. Adrian Korban, adrian3@windowslive.com		
// *************************************************************************************************
function cycgen(Rk,gg)
  n:=2;
  M:=RMatrixSpace(Rk,n,n)!0;
     for k:=1 to n do
       M[k]:=gg;
       temp:=gg;
       for t:=1 to (n-1) do
         temp[t+1]:=gg[t];
       end for;
       temp[1]:=gg[n];
       gg:=temp;
     end for; 
     return M;
end function;
function revcycgen(Rk,gg)
  n:=2;
  M:=RMatrixSpace(Rk,n,n)!0;
     for k:=1 to n do
       M[k]:=gg;
       temp:=gg;
       for t:=2 to n do
         temp[t-1]:=gg[t];
       end for;
       temp[n]:=gg[1];
       gg:=temp;
     end for; 
     return M;
  end function;
function GeneratorSigma(Rk,i)
	v1:=RSpace(Rk,2)![i[1],i[2]];
	v2:=RSpace(Rk,2)![i[3],i[4]];
	v3:=RSpace(Rk,2)![i[5],i[6]];
	v4:=RSpace(Rk,2)![i[7],i[8]];
	v5:=RSpace(Rk,2)![i[9],i[10]];
	v6:=RSpace(Rk,2)![i[11],i[12]];
	v7:=RSpace(Rk,2)![i[13],i[14]];
	v8:=RSpace(Rk,2)![i[15],i[16]];
	a1:=cycgen(Rk,v1);
	b1:=cycgen(Rk,v2);
	a2:=cycgen(Rk,v3);
	b2:=cycgen(Rk,v4);
	a3:=cycgen(Rk,v5);
	b3:=cycgen(Rk,v6);
	a4:=cycgen(Rk,v7);
	b4:=cycgen(Rk,v8);
	A1:=BlockMatrix(2,2,
	[
	a1,b1,
	Transpose(b1),Transpose(a1)
	]);
	A2:=BlockMatrix(2,2,
	[
	a2,b2,
	Transpose(b2),Transpose(a2)
	]);
	A3:=BlockMatrix(2,2,
	[
	a3,b3,
	Transpose(b3),Transpose(a3)
	]);
	A4:=BlockMatrix(2,2,
	[
	a4,b4,
	Transpose(b4),Transpose(a4)
	]);
	CMTemp:=BlockMatrix(4,4,
	[
	A1,A2,A3,A4,
	A2,A1,A4,A3,
	A3,A4,A1,A2,
	A4,A3,A2,A1
	]);
	return CMTemp;				
end function;  
FD:="~/Desktop/DIH16_DNA/ResultsDIH16_DNA.txt";
Codelist:={};
CodeDNA:=[];								
DDM:=0;
Rk<w>:=GF(4);
nVar:=16;
   Counter:=1;
   while true do
	xx:=RandomMatrix(Rk,1,nVar);
	CM:=GeneratorSigma(Rk,xx[1]);
	C:=LinearCode(CM);
	DM:=MinimumWeight(C);
	if (DM ge DDM) then
		CWD:=WeightDistribution(C);
		CWE<a,b,b,a>:=CompleteWeightEnumerator(C);
		if ([Length(C),Dimension(C),DM] notin Codelist) then
			printf "%o. Code-DNA is generated and save file \n", Counter;
			CodeDNA[Counter]:=C;
			T:=0;
			for i:=1 to #CWD do
				T:=T+CWD[i,2];
			end for;				
			AMG:=AutomorphismGroup(C);
			Write(FD,"*******************************************************");
			Write(FD,Counter);
			Write(FD,"*******************Minimum Weight**********************");
			Write(FD,DM);
			Write(FD,"*******************Dimension**********************");
			Write(FD,Dimension(C));
			Write(FD,"*******************Length**********************");
			Write(FD,Length(C));
			CodeParameters:=[Length(C),Dimension(C),DM];
			Write(FD,CodeParameters);										
			Write(FD,"***************Generator Matrix************************");
			Write(FD,CM);
			Write(FD,"-------------- Automorphism Group -----------------------");
			Write(FD,AMG);				
			Write(FD,"--------------- Weight Distribution --------------------");
			Write(FD,CWD);
			Write(FD,"--------------  GC Weight Enumerator --------------------");
			Write(FD,CWE);				
			Write(FD,"------------------ Number Of Words ----------------------");
			Write(FD,T);
			Codelist:=Codelist join {[Length(C),Dimension(C),DM]};				
			Counter:=Counter+1;
		end if;					
	end if;	
   end while;								
