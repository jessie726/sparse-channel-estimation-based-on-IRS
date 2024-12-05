function   kr = KhatriRao(F,G)     
		nR_F = size(F,1);        %F的行数
	  nR_G = size(G,1);      %G的行数
		mul = ones(nR_G,1);
		FF = kron(F,mul);     %通过kron函数实现对F矩阵的扩充，得到FF矩阵
		GG = repmat(G,nR_F,1);%通过repmat函数实现对G矩阵的扩充，得到GG矩阵
		kr = FF.*GG;
end