function   kr = KhatriRao(F,G)     
		nR_F = size(F,1);        %F������
	  nR_G = size(G,1);      %G������
		mul = ones(nR_G,1);
		FF = kron(F,mul);     %ͨ��kron����ʵ�ֶ�F��������䣬�õ�FF����
		GG = repmat(G,nR_F,1);%ͨ��repmat����ʵ�ֶ�G��������䣬�õ�GG����
		kr = FF.*GG;
end