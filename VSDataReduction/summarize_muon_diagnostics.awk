#!/usr/bin/awk -f
$10>10{
  for(i=0;i<4;i++){
    n[i]+=$(54+14*i);
    s[i]+=$(54+14*i)*$(55+14*i);
  }
};
END{
  for(i=0;i<4;i++){
    print n[i],s[i],s[i]/n[i];
  }
};
