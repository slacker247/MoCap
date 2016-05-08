fn SetObjectRotation obj rx ry rz =
	(
		-- Reset the object's transformation matrix so that
		-- it only includes position and scale information.
		-- Doing this clears out any previous object rotation.
		local translateMat = transMatrix obj.transform.pos
		local scaleMat = scaleMatrix obj.transform.scale
		obj.transform = scaleMat * translateMat

		-- Perform each axis rotation individually
		rotate obj (angleaxis rx [1,0,0])
		rotate obj (angleaxis ry [0,1,0])
		rotate obj (angleaxis rz [0,0,1])
	)


--in_file = openFile "Z:\\Software Dev\\SampleCode\\SerialPort\\Release\\Test2.txt"
in_file = openFile "Z:\\Software Dev\\Embedded\\MoCap\\Test.txt.ft"

ax = #()
ay = #()
az = #()
timeIndex = #()

index = 1

while not eof in_file do
(
	lne = readline in_file

	if lne == #eof then
		break

	vIndex = 0
	test = true
	while test do
	(
		pos = findString lne " "

		if pos != undefined then
		(
			val = substring lne 1 (pos - 1)
			if pos == 1 then
			(
				lne = substring lne (pos+1) lne.count
				continue
			)
			
			if vIndex == 0 then
				append timeIndex val
			if vIndex == 1 then
				append ax val
			if vIndex == 2 then
				append ay val
			if vIndex == 3 then
				append az val
			lne = substring lne (pos+1) lne.count
			vIndex = vIndex + 1
		)
		else
		(
			append az lne
			test = false
		)
	)
	--format "%, %, %, %\n" timeIndex[index] ax[index] ay[index] az[index]
	index = index + 1
)

close in_file

s = sphere()
s.radius = 10.0

translateMat = transMatrix s.transform.pos
scaleMat = scaleMatrix s.transform.scale
s.transform = scaleMat * translateMat

index = 1
x = ax[index]
y = ay[index]
z = az[index]

animationRange = interval 0f (timeIndex[timeIndex.count] as time)

at time ((timeIndex[index] as time) - (timeIndex[1] as time)) animate on s.rotation.x_rotation = (ax[index] as Number)
at time ((timeIndex[index] as time) - (timeIndex[1] as time)) animate on s.rotation.y_rotation = (ay[index] as Number)
at time ((timeIndex[index] as time) - (timeIndex[1] as time)) animate on s.rotation.z_rotation = (az[index] as Number)
lastTime = ("0:0.0" as time)

for i = 0 to ax.count do
(
	if ax[index] == undefined  OR
	   ay[index] == undefined  OR
	   az[index] == undefined  then
		continue

/*	at time ((timeIndex[index] as time) - (timeIndex[1] as time)) animate on SetObjectRotation s (ax[index] as Float) (ay[index] as Float) (az[index] as Float) 
--	s.transform = scaleMat * translateMat
	s.rotation = (eulerangles 0 0 0)*/
		at time ((timeIndex[index] as time) - (timeIndex[1] as time)) animate on s.rotation = ( eulerangles ((ax[index] as Float)*3.14/180) ((ay[index] as Float)*3.14/180) ((az[index] as Float)*3.14/180) )
/*
	if lastTime == undefined OR
		(lastTime + 10f) < (((timeIndex[index] as time) - (timeIndex[1] as time)) - lastTime)
	    then
	(
		s.rotation.x_rotation = 0;
		s.rotation.y_rotation = 0;
		s.rotation.z_rotation = 0;
		
		at time ((timeIndex[index] as time) - (timeIndex[1] as time)) animate on s.rotation.x_rotation = (ax[index] as Number)
		at time ((timeIndex[index] as time) - (timeIndex[1] as time)) animate on s.rotation.y_rotation = (ay[index] as Number)
		at time ((timeIndex[index] as time) - (timeIndex[1] as time)) animate on s.rotation.z_rotation = (az[index] as Number)
		if lastTime == undefined then
		(
			lastTime = ((timeIndex[index] as time) - (timeIndex[1] as time))
		)
		else
		(
			lastTime = ((timeIndex[index] as time) - (timeIndex[1] as time)) - lastTime
		)
	)*/
	index = index + 1
)
