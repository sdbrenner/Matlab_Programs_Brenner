function jd=julday(x)
%converts date to julian date; 0000 on jan 1=1.000
cov1=[0 31 59 90 120 151 182 213 244 274 305 335 ];
cov2=[0 31 60 91 121 152 183 214 245 275 306 336 ];
yr=x(1);
mon=x(2);
day=x(3);
hr=x(4);
min=x(5);
sec=x(6);
if fix(yr/4)==0
	jd=cov2(mon)+day+hr/24+min/24/60+sec/24/3600;
else
	jd=cov1(mon)+day+hr/24+min/24/60+sec/24/3600;
end;