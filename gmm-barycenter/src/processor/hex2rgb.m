function rgb = hex2rgb(hexString)
	if size(hexString,2) ~= 7
		error('invalid input: not 7 characters');
	else
		r = double(hex2dec(hexString(2:3)))/255;
		g = double(hex2dec(hexString(4:5)))/255;
		b = double(hex2dec(hexString(6:7)))/255;
		rgb = [r, g, b];
	end
end