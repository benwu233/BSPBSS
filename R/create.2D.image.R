create.circle.in.2D.image = function(voxels,center = c(0.5,0.5),radius = 0.2){
  return((voxels[,1]-center[1])^2+(voxels[,2]-center[2])^2<radius^2)
}

create.triangle.in.2D.image = function(voxels,x0,y0, x1,y1, x2,y2){
  return(points.in.triangle(voxels[,1],voxels[,2],x0,y0, x1,y1, x2,y2))
}

create.square.in.2D.image = function(voxels,x0,y0,x1,y1){
  return(points.in.rectangle(voxels[,1],voxels[,2],x0,y0,x1,y1))
}


points.in.rectangle = function(x,y,x0,y0,x1,y1){
  return((x>x0 & x<x1) & (y>y0 & y<y1))
}

points.in.triangle = function(x,y, x0,y0, x1,y1, x2,y2) {
  s = y0 * x2 - x0 * y2 + (y2 - y0) * x + (x0 - x2) * y;
  t = x0 * y1 - y0 * x1 + (y0 - y1) * x + (x1 - x0) * y;

  res = rep(FALSE,length=length(x))

  false_idx = which((s < 0) != (t < 0))
  if(length(false_idx)<length(x)){
    A = -y1 * x2 + y0 * (x2 - x1) + x0 * (y1 - y2) + x1 * y2;
    if (A < 0.0) {
      s = -s;
      t = -t;
      A = -A;
    }
    res = (s > 0 & t > 0 & (s + t) <= A)
  }
  res[false_idx] = FALSE
  return(res)
}
