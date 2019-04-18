SetDirectory["$(pwd)"]

mx=mus[[1,1]];
my=mus[[1,2]];

day:=DateList[][[3]];
month:=DateList[][[2]];
hour:=DateList[][[4]];
minute:=DateList[][[5]];

Table[
	Module[
		{rev,ref,ndestats},
		ToString[StringForm["Starting NDE polar calculation, \[Epsilon] = `5` ; `1`/`2` at `3`h`4`m",month,day,hour,minute,N[eps]]]>>>mylog.txt;
		{{rev,ref},ndestats}=med[Unevaluated[ndespolar[10,mx,my,VKeld,4.89,chiphos,0,eps[[1]],10^12,None]]];
		ToString[StringForm["NDE calculation finished, `1`/`2` at `3`h`4`m",month,day,hour,minute]]>>>mylog.txt;
		Export[ToString@StringForm["ndepolar.m",eps[[2]]],{rev,ref}];
		Export[ToString@StringForm["ndestats.m",eps[[2]]],ndestats];
		ToString[StringForm["Data saved, looping back... `1`h`2`m",hour,minute]]>>>mylog.txt;
	],
	{eps,{{10^-3,"1e3"},{9*10^-4,"9e4"},{8*10^-4,"8e4"},{7*10^-4,"7e4"},{6*10^-4,"6e4"}}}
]
