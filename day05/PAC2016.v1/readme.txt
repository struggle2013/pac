Э�����㷨�Ż�

1. ������ܣ�
�����ṩ4��Դ�ļ���main.f90, check.f90, Covariance.f90��Covariance1.f90���Լ�һ��Makefile�ļ���
���У�

a) main.f90�ṩ���ݳ�ʼ�������ù��ܼ��㺯���������֤��

b) check.f90����У�顣��check.f90�ļ�IsMissingPheno�����е��ж���������ı䡣
c) Covariance1.f90��Covariance.f90���ǹ��ܼ��㺯����Ŀǰ�����߼�ʵ��һ�£�������֤ʹ�á�


ֱ��make�������ɿ�ִ���ļ� main��ֱ�����п�ִ���ļ��ó������ʱ��
���룺
$ make
ifort -c main.f90
ifort -c Covariance1.f90
ifort -c Covariance.f90
ifort -c check.f90
ifort  main.o Covariance1.o Covariance.o check.o -o main

���У�
./main
 total time(s):    ����ʱ��    
 verification: correct

2. ��������
a) �����޸� main.f90, check.f90 �� Covariance1.f90 �����ļ������޸ģ������ڳɼ�ֱ�����ϡ�
b) �����Գ����ڵ���Covariance����ʱ����������ֵֻ����һ�ֳ��������
ѡ����Ҫ���ݳ����߼������Ż�������ֻ��Դ��ݵĲ���ֵ���ض��Ż������򱾻��ڳɼ����ϡ�
c) ��������֤ͨ�����ɼ���Ч�����򱾻��ڳɼ����ϡ�
d) �����ṩ����ƽ̨KNL������ʹ��ƽ̨Ԥ��װ��Intel compiler 2016.0.0.3.
e) �����Ż��Ƚϣ��Գ��������ʱΪ�������ݣ���ʱԽ���߳ɼ�Խ�ߡ�


