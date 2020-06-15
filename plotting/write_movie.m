function write_movie(F, file_name, frame_rate)

writerObj = VideoWriter(file_name,'MPEG-4');
writerObj.FrameRate = frame_rate;

open(writerObj);
writeVideo(writerObj, F);
close(writerObj);

end