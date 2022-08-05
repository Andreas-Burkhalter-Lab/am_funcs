function result=smooth(data,windowsize)
b=ones(1,windowsize)/windowsize;
a=1;
result=filter(b,a,data);
