function vec_mrg = tdrVectorMerge(vec1,vec2)


vec_mrg = vec1;
vec_mrg.name = [vec1.name;vec2.name];

nt1 = size(vec1.response,2);
nt2 = size(vec2.response,2);

if nt1==1 && nt2==1
    vec_mrg.response = cat(3,vec1.response,vec2.response);
elseif nt1>1 && nt2==1
    rr1 = vec1.response;
    rr2 = repmat(vec2.response,[1 nt1 1]);
    vec_mrg.response = cat(3,rr1,rr2);
elseif nt1==1 && nt2>1
    rr1 = repmat(vec1.response,[1 nt2 1]);
    rr2 = vec2.response;
    vec_mrg.response = cat(3,rr1,rr2);
else
    error('wrong vector size');
end
    

