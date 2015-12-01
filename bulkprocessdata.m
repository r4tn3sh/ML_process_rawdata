clc;
clear;
close all;

% PARAMETERS NEEDED TO BE ADJUSTED
for ia = 10:10
    for ib = 4:4:20
        for ic = ia:320
            if rem(ia*ib+ic,80)==0
                foldername = strcat('/home/ratensh/ONRDATACOLLECTION/datacollection/1_pair/sync_data/',num2str(ia),'_',num2str(ib),'_',num2str(ic),'/');
                command = ['cp ',foldername,'* ','/home/ratensh/simulationORBITdata/data/'];
                [status,cmdout] = system(command);
                command = '/home/ratensh/simulationORBITdata/data/generatemat';
                [status,cmdout] = system(command);
                configdata=textread('/home/ratensh/simulationORBITdata/data/common_config','%s');
                configlen = length(configdata);
                for i=1:configlen
                    if strcmp(configdata{i},'mode')==1
                        mode=str2num(configdata{i+1});
                    elseif strcmp(configdata{i},'shseq_len')==1
                        shseq_len=str2num(configdata{i+1});
                    elseif strcmp(configdata{i},'shseq_rep')==1
                        shseq_rep=str2num(configdata{i+1});
                    elseif strcmp(configdata{i},'lgseq_len')==1
                        lgseq_len=str2num(configdata{i+1});
                    elseif strcmp(configdata{i},'dataseq_len')==1
                        dataseq_len=str2num(configdata{i+1});
                    end
                end
                seq_len = shseq_len*shseq_rep+lgseq_len;
                totallen=1280 + seq_len;
                load('/home/ratensh/simulationORBITdata/data/longdata.mat')
                load('/home/ratensh/simulationORBITdata/data/lgthdata.mat')
                load('/home/ratensh/simulationORBITdata/data/corrdata.mat')
                load('/home/ratensh/simulationORBITdata/data/gaindata.mat')
                load('/home/ratensh/simulationORBITdata/data/orgdata.mat')
                configdata=textread('/home/ratensh/simulationORBITdata/data/common_config','%s');
                configlen = length(configdata);
                for i=1:configlen
                    if strcmp(configdata{i},'mode')==1
                        mode=str2num(configdata{i+1});
                        break;
                    end
                end


                datalen = length(corrdata);
                dataoffset = ceil(0.1*datalen); 
                totalpeaks = sum(corrdata(dataoffset:datalen));
                previndex=0;
                currindex=0;
                initfalsecnt=0;
                % Find the location of potention correct peaks and discard the wrong
                for i=totallen+dataoffset:datalen
                    pflag(i)=0;
                    if corrdata(i)==1
                        for j=i-totallen:-totallen:1
                            if corrdata(j)==1
                                pflag(i)=1;
                                break;
                            end
                        end
                        if pflag(i)==0
                            initfalsecnt = initfalsecnt+1;
                        end
                    end
                end
                doublepeakcnt=0;
                missedpeakcnt=0;
                for i=dataoffset:datalen
                    indexdiff(i)=0;
                    flag(i)=0;
                    flagval(i)=0;
                    correctpeakflag(i) = 0;
                    if pflag(i)==1
                        currindex = i;
                        if previndex>0
                            indexdiff(i) = currindex-previndex;
                            if rem(indexdiff(i),totallen)>1 && rem(indexdiff(i),totallen)<totallen-1
                                flag(i)=1;
                                flagval(i)=indexdiff(i);
                            elseif i>1 && rem(indexdiff(i-1),totallen)==totallen-1 && rem(indexdiff(i),totallen)~=1
                                flag(i)=1;
                                flagval(i)=indexdiff(i);
                            elseif i>1 && rem(indexdiff(i-1),totallen)~=totallen-1 &&...
                                    rem(indexdiff(i-1),totallen)>0 && rem(indexdiff(i),totallen)==1
                                flag(i)=1;
                                flagval(i)=indexdiff(i);
                            elseif i>1 && rem(indexdiff(i-1),totallen)>0 && rem(indexdiff(i),totallen)==1
                                correctpeakflag(i) = 1;
                                doublepeakcnt=doublepeakcnt+1;
                                missedpeakcnt = missedpeakcnt+(indexdiff(i-1)/totallen)-1;
                            elseif rem(indexdiff(i),totallen)==0
                                correctpeakflag(i) = 1;
                                missedpeakcnt = missedpeakcnt+(indexdiff(i)/totallen)-1;
                            else
                                a=0;
                            end
                %             if i>1 && rem(indexdiff(i),totallen)==1 &&
                %                 flag(i)=1;
                %                 flagval(i)=indexdiff(i);
                %             end
                        end
                    end
                    previndex=currindex;
                end
                falsepeakcnt=sum(flag)+initfalsecnt;
                correctpeakcnt = sum(correctpeakflag);
                maxcorr = 1.2*max(abs(longdata));
                figure;
                hold on;
                %plot(maxcorr*corrdata,'g');
                plot(maxcorr*flag,'r');
                plot(maxcorr*correctpeakflag,'g');
                plot(abs(longdata),'b');
                plot(lgthdata,'c');
                title(['L_{sh}=',num2str(ia),', K=',num2str(ib),', L_{lg}=',num2str(ic)])

                str = sprintf('*****************%d:%d:%d*******************',ia,ib,ic);
                disp(str);
                str = sprintf('False alarm rate: %f',falsepeakcnt/(totalpeaks-doublepeakcnt));
                disp(str);
                str = sprintf('Missed rate: %f',missedpeakcnt/((datalen-dataoffset)/totallen));
                disp(str);
            end
        end
    end
end