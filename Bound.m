
function s = Bound(s, Lb, Ub,best1,best2) 
temp = s;
III = temp < Lb;
temp(III)=best2(III)+rand*(best1(III)-best2(III));
J = temp > Ub;
temp(J)=best2(J)+rand*(best1(J)-best2(J));
s = temp;
end