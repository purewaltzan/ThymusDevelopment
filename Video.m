%%loading files list and select
files=ls;
k=zeros(numel(files),1)==1;
%test whether the files are merged image or specific channel
for i=1:size(files,1)
    k(i)=contains(files(i,:),'.tif ');
end
selected_files=files(k,:);
%$generate video
outputVideo = VideoWriter(fullfile('Video.avi'));
outputVideo.FrameRate = 4;
open(outputVideo)
for ii = 1:sum(k)
    img = imread(fullfile(selected_files(ii,:)));
    writeVideo(outputVideo,img);
end
close(outputVideo)