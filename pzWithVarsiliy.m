function [P,Z,dist] = pzWithVarsiliy(H,rate)
r = 0;
tmp = 0;
P = [];
Z = [];
while r<rate
  r = r + H(tmp+256);
  P = [P,tmp];
  if tmp < 0
      tmp =  -tmp;
  else
      tmp = -tmp - 1;
  end
end
P = sort(P);
Tn = min(P);
Tp = max(P);
Z = [];
tmp = 1;
while H(tmp) == 0
    tmp=tmp+1;
end
tmp = tmp - 1;
if Tn ~= 0
    Z = tmp+Tn+1:1:tmp;
    Z = Z - 256;
end
tmp = 511;
while H(tmp) == 0
    tmp=tmp-1;
end
tmp = tmp+1;
tz = tmp:1:tmp+Tp;
tz = tz-256;
Z = [Z,tz];

[Rate,dist] = rateAndDist(P+256,Z+256,H);

