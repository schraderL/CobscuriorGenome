BEGIN{ FS=OFS="," }
NR==1{
   header=$0
   n=split("Dup,Counter",h)
   for (i=1; i<=NF; i++)
      for (j=1; j<=n; j++) header=header OFS $i"_"h[j]
   printf("%s,EntireLine_Dup,EntireLine_Counter\n", header)
   next
}
{
    r[++lines]=$0
    for (col=1; col<=NF; col++) v[col][$col]++
    v[col][$0]++
}
# END {
#    for (l=1; l<=lines; l++){
#       n=split(r[l], s)
#       res=""
#       for (c=1; c<=n; c++)
#          res=res OFS output(v,c,s[c])
#       res=res OFS output(v,c,r[l])
#       print r[l] res
#    }
# }
# function output(arr, col, val){
#     return sprintf("%s,%s", (arr[col][val] > 1? "Yes" : "No"), ++count[col][val])
# }
