export NP2SP
function NP2SP(hr::HR, orbital::ORBITAL;  buildindex = 'Y')

	spinhr = spinHR(hr; buildindex)
	spinorbital = spinORBITAL(orbital)

	return spinhr, spinorbital
end


# function NP2SP(hr::HR, orbital::ORBITAL;  buildindex = 'Y')

# 	spinhr = spinHR(hr; buildindex)
# 	spinorbital = spinORBITAL(orbital)

# 	return spinhr, spinorbital
# end