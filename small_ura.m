function sURA=small_ura(xapp,yapp)

x= -abs(xapp):abs(xapp);
y= -abs(yapp):abs(yapp);
[xx, yy]= meshgrid(x,y);
xx2=reshape(xx,[],1);
yy2=reshape(yy,[],1);
sURA= [xx2,yy2];







end 