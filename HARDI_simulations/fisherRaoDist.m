function y=fisherRaoDist (odf1, odf2)
sqOdf1 = sqrt(odf1/sum(odf1));
sqOdf2 = sqrt(odf2/sum(odf2));
y = acos(dot(sqOdf1, sqOdf2));
end