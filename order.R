setwd("/home/hou/test/1111")


files<-list.files("/home/hou/test/1111")
par<-"wrong"
for (i in 1:length(files))
{
if (files[i]=="ref.fasta")
{par<-"right";break}
i<-i+1}
if (par=="right")
{
seq_line<-1
ref<-read.csv("ref.fasta",header=F)
ref<-as.matrix(ref)
ref_char<-unlist(strsplit(ref,split=""))
ref_length<-length(ref_char)

#互补
complement<-function(seq)
{
if (seq=="A") {seq1<-"T"}
if (seq=="a") {seq1<-"t"}
if (seq=="C") {seq1<-"G"}
if (seq=="c") {seq1<-"g"}
if (seq=="G") {seq1<-"C"}
if (seq=="g") {seq1<-"c"}
if (seq=="T") {seq1<-"A"}
if (seq=="t") {seq1<-"a"}
if (seq=="|") {seq1<-"|"}
if (seq=="-") {seq1<-"-"}
if (seq==" ") {seq1<-" "}
if (seq=="") {seq1<-""}
seq1
}

#去掉没有匹配上的文件
aln_name<-matrix(nrow=length(files),ncol=1)
name<-1;num_normal<-1
for (i in 1:length(files))
{
len<-nchar(files[i])
if (substr(files[i],len-3,len)==".aln")
{
csv<-read.csv(files[i],header=T,sep=",",skip=14)
csv<-as.matrix(csv)
for (m in 1:nrow(csv))
{
if (csv[m,1]=="***** No hits found *****")
{break}
else {m<-m+1}}
if (m==nrow(csv)+1) {aln_name[name,]<-files[i];num_normal<-num_normal+1;name<-name+1}
i<-i+1}
else {i<-i+1}}
if (num_normal>1)
{
#对匹配的文件进行排序
if (name==2)
{aln_name<-rbind(aln_name[1],aln_name[1]);name<-3}
order_ref<-matrix(nrow=name-1,ncol=3)
m<-1

for (line in 1:(name-1))
{
csv<-read.csv(aln_name[line],header=T,sep=",",skip=14)
csv<-as.matrix(csv)
for (i in 1:nrow(csv))
{
if (substr(csv[i,1],1,13)=="Lambda      K")
{
end_line<-i-1;break}
else {i<-i+1}}
for (i in end_line:end_line)
{
if (substr(csv[i,1],1,5)=="Sbjct")
{
for (j in nchar(csv[i,1]):(nchar(csv[i,1])-10))
{
if ((substr(csv[i,1],j,j)=="A") || (substr(csv[i,1],j,j)=="C") || (substr(csv[i,1],j,j)=="G") || (substr(csv[i,1],j,j)=="T"))
{Sbjct_end<-as.numeric(substr(csv[i,1],j+1,nchar(csv[i,])));order_ref[m,3]<-Sbjct_end;order_ref[m,1]<-aln_name[line];break}
else {j<-j-1}};break}
i<-i-1}
for (i in 1:nrow(csv))
{
if (substr(csv[i,1],1,5)=="Query") {Sbjct_start<-as.numeric(substr(csv[i+2,1],6,12));order_ref[m,2]<-Sbjct_start;break}
else {i<-i+1}}
line<-line+1;m<-m+1}

#排序
for (i in 1:nrow(order_ref))
{
start<-as.numeric(order_ref[i,2])
end<-as.numeric(order_ref[i,3])
if (start>end) {x<-start;y<-end;order_ref[i,2]<-end;order_ref[i,3]<-start}
i<-i+1}
aln_name_order<-order_ref[order(as.numeric(order_ref[,2])),]

final<-matrix(nrow=2*name,ncol=2)
#参考序列的行生成
query_kong<-""
file_length<-3
for (m in 1:(60-file_length))
{
query_kong<-paste(query_kong," ",sep="")
m<-m+1}
final[1,1]<-paste(query_kong,"Num ",sep="")
final[2,1]<-paste(query_kong,"Ref ",sep="")
num<-"";x<-""
for(i in 1:ref_length) {final[2,2]<-paste(x,ref_char[i],sep="");x<-final[2,2];i<-i+1}

number<-""
for (i in 1:nchar(final[2,2])) 
{
if ((i%%10==0) && (i>=10) && (i<100))
{
final[1,2]<-paste(number,i,"        ",sep="");number<-final[1,2];i<-i+1
}
if ((i%%10==0) && (i>=100) && (i<1000))
{
final[1,2]<-paste(number,i,"       ",sep="");number<-final[1,2];i<-i+1
}
if ((i%%10==0) && (i>=1000) && (i<10000))
{
final[1,2]<-paste(number,i,"      ",sep="");number<-final[1,2];i<-i+1
}
if ((i%%10==0) && (i>=10000) && (i<100000))
{
final[1,2]<-paste(number,i,"     ",sep="");number<-final[1,2];i<-i+1
}
i<-i+1
}
x<-paste("1","        ",final[1,2],sep="")
final[1,2]<-x
for (i in 1:nrow(csv))
{
if (substr(csv[i,1],1,5)=="Query") {sj_st<-as.numeric(substr(csv[i+2,1],6,12));sj_st_line<-i;break}
else {i<-i+1}}

#匹配的行生成
k<-0


output_insert<-matrix(nrow=nrow(aln_name_order)+1,ncol=800)
output_del<-matrix(nrow=nrow(aln_name_order)+1,ncol=800)
output_mut<-matrix(nrow=nrow(aln_name_order)+1,ncol=800)




for (line in 1:nrow(aln_name_order))
{
output_del[line+1,1]<-aln_name_order[line,1];output_mut[line+1,1]<-aln_name_order[line,1];output_insert[line+1,1]<-aln_name_order[line,1]
csv<-read.csv(aln_name_order[line,1],header=T,sep=",",skip=14)
csv<-as.matrix(csv)
#先找到文件的头和尾
#Query和ref的起始位点，开始行
for (i in 1:nrow(csv))
{
if (substr(csv[i,1],1,5)=="Query")
{
for (j in 1:30)
{
if ((substr(csv[i,1],j,j)=="A") || (substr(csv[i,1],j,j)=="C") || (substr(csv[i,1],j,j)=="G") || (substr(csv[i,1],j,j)=="T"))
{
Q_st<-as.numeric(substr(csv[i,1],6,j-1));R_st<-as.numeric(substr(csv[i+2,1],6,j-1));Q_st_site<-j;break}
j<-j+1};cut_st<-i;break}
else {i<-i+1}}

for (i in cut_st:nrow(csv))
{
if (substr(csv[i,1],1,6)==" Score")
{cut_end<-i-1;break}
if (substr(csv[i,1],1,13)=="Lambda      K")
{cut_end<-i-1;break}
else {i<-i+1}
}

#把序列赋值给一个新的矩阵
new<-csv[cut_st:cut_end,]
end<-length(new)


#Query和ref的结束位点吗，结束行
for (j in nchar(new[end]):1)
{
if ((substr(new[end],j,j)=="A") || (substr(new[end],j,j)=="C") || (substr(new[end],j,j)=="G") || (substr(new[end],j,j)=="T"))
{R_end_site<-j;break}
else {j<-j+1}}

Q_end<-as.numeric(substr(new[end-2],R_end_site+1,nchar(new[end-2])))
R_end<-as.numeric(substr(new[end],R_end_site+1,nchar(new[end])))



#前两行生成
query_kong<-""
file_length<-nchar(aln_name_order[line,1])
for (m in 1:(60-file_length))
{
query_kong<-paste(query_kong," ",sep="")
m<-m+1}
query_kong1<-""
for (m in 1:61)
{
query_kong1<-paste(query_kong1," ",sep="")
m<-m+1}

final[2*(line+1),1]<-paste(query_kong,aln_name_order[line,1]," ",sep="")

lable<-paste("Ref_",Q_st,":",Q_end,"/Seq_",R_st,":",R_end,sep="")
for (m in 1:(60-(nchar(lable))))
{
lable<-paste(" ",lable,sep="")
m<-m+1}
final[2*(line+1)-1,1]<-paste(lable," ",sep="")

#合并文件为单行
merge<-matrix(nrow=3,ncol=1)
merge[1,1]<-"";merge[2,1]<-"";merge[3,1]<-""
for (i in 1:(length(new)/3-1))
{
line1<-substr(new[3*i-2],Q_st_site,Q_st_site+59)
line2<-substr(new[3*i-1],Q_st_site,Q_st_site+59)
line3<-substr(new[3*i],Q_st_site,Q_st_site+59)
merge[1,1]<-paste(merge[1,1],line1,sep="")
merge[2,1]<-paste(merge[2,1],line2,sep="")
merge[3,1]<-paste(merge[3,1],line3,sep="")
i<-i+1}
merge[1,1]<-paste(merge[1,1],substr(new[3*i-2],Q_st_site,R_end_site),sep="")
merge[2,1]<-paste(merge[2,1],substr(new[3*i-1],Q_st_site,R_end_site),sep="")
merge[3,1]<-paste(merge[3,1],substr(new[3*i],Q_st_site,R_end_site),sep="")





#判断ref的方向,query的方向
if (R_end > R_st)
{
insert<-1;del<-1;mut<-1
seq_len<-nchar(merge[3,1])
query_kong1<-"";query_kong2<-""
if (R_st==1) {query_kong1<-"";query_kong2<-""}
else
{
for (m in 1:(R_st-1))
{
query_kong1<-paste(query_kong1," ",sep="")
query_kong2<-paste(query_kong2," ",sep="")
m<-m+1}
}
final[2*(line+1)-1,2]<-query_kong1
final[2*(line+1),2]<-query_kong2


#找到最终序列的del，insert等
for (j in 1:seq_len)
{
if (substr(merge[3,1],j,j)!="-")
{
final[2*(line+1)-1,2]<-paste(query_kong1,substr(merge[2,1],j,j),sep="");query_kong1<-final[2*(line+1)-1,2]
final[2*(line+1),2]<-paste(query_kong2,substr(merge[1,1],j,j),sep="");query_kong2<-final[2*(line+1),2]
#del_information
output_del[1,1]<-"Del_information"
if (substr(merge[1,1],j,j)=="-") 
{
Del_Q<-Q_st+j-1;Del_R<-R_st+j-1;
output_del[line+1,del+1]<-paste(Del_Q,":",Del_R);del<-del+1}

#mut_information
output_mut[1,1]<-"Mut_information"
if ((substr(merge[1,1],j,j)!=substr(merge[3,1],j,j)) && (substr(merge[1,1],j,j)!="-"))
{
M_Q<-Q_st+j-1;M_R<-R_st+j-1;M_BASE_Q<-substr(merge[1,1],j,j);
M_BASE_R<-substr(merge[3,1],j,j)
output_mut[line+1,mut+1]<-paste(M_Q,":",M_R,"_",M_BASE_Q,":",M_BASE_R);mut<-mut+1
}
}
else 
{
#insert_information
output_insert[1,1]<-"Insert_information"
I_Q<-Q_st+j-1;I_R<-R_st+j-1
output_insert[line+1,insert+1]<-paste(I_Q,":",I_R);insert<-insert+1
}
j<-j+1
}
}


#ref是反向的
else
{
insert<-1;del<-1;mut<-1
seq_len<-nchar(merge[3,1])
#标准ref生成
query_kong1<-"";query_kong2<-""
if (R_end==1) {query_kong1<-"";query_kong2<-""}
else
{
for (m in 1:(R_end-1))
{
query_kong1<-paste(query_kong1," ",sep="")
query_kong2<-paste(query_kong2," ",sep="")
m<-m+1}
}
final[2*(line+1)-1,2]<-query_kong1
final[2*(line+1),2]<-query_kong2

#找到最终序列的del，insert等
for (j in 1:seq_len)
{
if (substr(merge[3,1],seq_len-j+1,seq_len-j+1)!="-")
{
final[2*(line+1)-1,2]<-paste(query_kong1,complement(substr(merge[2,1],seq_len-j+1,seq_len-j+1)),sep="");query_kong1<-final[2*(line+1)-1,2]
final[2*(line+1),2]<-paste(query_kong2,complement(substr(merge[1,1],seq_len-j+1,seq_len-j+1)),sep="");query_kong2<-final[2*(line+1),2]
#del_information
output_del[1,1]<-"Del_information"
if (substr(merge[1,1],seq_len-j+1,seq_len-j+1)=="-") 
{
Del_Q<-Q_end-j+1;Del_R<-R_end+j-1;
output_del[line+1,del+1]<-paste(Del_Q,":",Del_R);del<-del+1}

#mut_information
output_mut[1,1]<-"Mut_information"
if ((substr(merge[1,1],seq_len-j+1,seq_len-j+1)!=substr(merge[3,1],seq_len-j+1,seq_len-j+1)) && (substr(merge[1,1],seq_len-j+1,seq_len-j+1)!="-"))
{
M_Q<-Q_end-j+1;M_R<-R_end+j-1;;M_BASE_Q<-substr(merge[1,1],seq_len-j+1,seq_len-j+1);
M_BASE_R<-complement(substr(merge[3,1],seq_len-j+1,seq_len-j+1))
output_mut[line+1,mut+1]<-paste(M_Q,":",M_R,"_",M_BASE_Q,":",M_BASE_R);mut<-mut+1
}
}
else 
{
#insert_information
output_insert[1,1]<-"Insert_information"
I_Q<-Q_end-j+1;I_R<-R_end+j-1;
output_insert[line+1,insert+1]<-paste(I_Q,":",I_R);insert<-insert+1
}
j<-j+1
}
}
line<-line+1}
write.table(final,"seq.txt",quote=FALSE,col.names=FALSE,row.names=FALSE,na=" ")
output<-rbind(output_del,output_insert,output_mut)
write.table(output,"site.txt",quote=FALSE,col.names=FALSE,row.names=FALSE,na=" ",sep="\t")
}
}

#生成pdf图片
#source("https://bioconductor.org/biocLite.R")
#biocLite("sangerseqR")
#library(sangerseqR)
#setwd("F://production_test/")
#files<-list.files("F://production_test/")
#for (m in 1:length(files))
#{
#if (substr(files[m],nchar(files[m])-2,nchar(files[m]))=="ab1")
#{
#jpeg(file=paste(files[m],".jpeg",sep=""))
#hetsangerseq <- readsangerseq(files[m])
#x<-chromatogram(hetsangerseq, width=100, height=2,showcalls='both', filename=paste(files[m],".pdf",sep=""))
#x<-chromatogram(hetsangerseq, width=100, height=2,showcalls='both')
#dev.off()
#m<-m+1}
#else {m<-m+1}}


