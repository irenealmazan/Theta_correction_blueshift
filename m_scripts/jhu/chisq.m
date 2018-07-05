function chi = chisq(array1, array2, sigma, normflag);

if(normflag)
    array1 = array1/max(array1);
    array2 = array2/max(array2);
end

array3 = real(array1) - real(array2);
array3=array3 / sigma;
array3= array3.^2;
chi=sum(array3);


