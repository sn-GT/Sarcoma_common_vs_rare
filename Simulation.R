library(ggplot2)
library(cowplot)
library(Rmisc)

n=100000000
tot <- n
prev <- 0.00007
smallprev <-0.000065
carrier_perc <- 0.02
allids <- data.frame("ID" = c(1:n), "RandNo" = runif(tot, min=0, max=1))

#10% carriers
allids$carriers <- 0
#allids$carriers[allids$RandNo <= carrier_perc] <- 1
# randomly assign
carriercount <- n*carrier_perc 
allids$carriers[sample(nrow(allids), carriercount)] <- 1

print(table(allids$carriers)[2]/tot)

#car_case_perc <- 0.01 # penetrance
clist <- c(0.0003, 0.00025, 0.0002, 0.00015, 0.0001)
fullres <- data.frame()

for(it in 1:5){
	print(it)
	resdf <- data.frame()
	for(car_case_perc in clist){
		print(car_case_perc)
		tot_cases <- prev*tot 
		carr_cases <- tot*carrier_perc*car_case_perc
		noncarr_cases <- tot_cases - carr_cases
		
		carriers <- allids[allids$carriers == 1,]
		c1 <- carriers[sample(nrow(carriers), carr_cases),]
		c1$CASE <- "CASE"
		c2 <- carriers[!carriers$ID %in% c1$ID,]
		c2$CASE <- "CTRL"
		
		c <- rbind(c1, c2)  
		nc <- allids[allids$carriers ==0,]
		nc$CASE <- "TBD"
		
		cb <- rbind(c, nc)
		
		# generate PRS norm dist
		cb$PRS <-  rnorm(runif(nrow(cb)),5, 0.75)
		# plot PRS distribution
		ggplot(cb, aes(x= PRS)) + geom_density() + theme_bw() 
		
		# get carrier cases and their prev in each PRS percentile bin 
		# prev - perc based on carriers cases
		n=100
		dx <- cb #[cb$carriers == 1,]
		names(dx)[5] <- "invGRS"
		df <- data.frame()
		dfull <- data.frame()
		for(i in 0:(n-1)){
			print(i)
			if(i == 0){
				dx1 <- dx[which(dx$invGRS >= quantile(dx$invGRS, i/n) & dx$invGRS <= quantile(dx$invGRS, (i+1)/n)),]
			}
			else {
				dx1 <- dx[which(dx$invGRS > quantile(dx$invGRS, i/n) & dx$invGRS <= quantile(dx$invGRS, (i+1)/n)),]
			}
			prev3 <- table(dx1$CASE)["CASE"]/nrow(dx1)
			if(is.na(prev3)){prev3 = 0}
			
			ncarrier_case <- table(dx1$CASE)["CASE"]
			if(is.na(ncarrier_case)){ncarrier_case = 0}
			
			# full dx 
			dx1$PGSbin <- i+1
			dx1$CarrCase <- ncarrier_case
			dx1$Prev <- prev3
			dx1$N <- nrow(dx1)
			dfull <- rbind(dfull, dx1)
			
			# prev perc df
			df2 = data.frame("PGS" = i+1, "nCarrierCase" =  ncarrier_case,"meanPGS" = mean(dx1$invGRS) , "Prev3" = prev3) 
			df <- rbind(df, df2)
		}
		p_carriers <- ggplot(df, aes(x=PGS, y=Prev3)) + 
		geom_point( size = 2.3) + theme_bw() + xlab("Percentile PRS") + ylab(paste0("Prevalence of disease")) + 
		geom_smooth(method = "lm",formula = y ~ poly(x, 3, raw = T),  size = 1, se = T) + #ylim(0.001, 0.0082) + 
		ggtitle("Carriers")
		p_carriers
		print("carrierdone")
		#################################################################
		# Polygenic background in full data
		
		# prev per bin by PGS background
		# 4x bottom vs top decile
		#p1 <- prev/2
		#p50 <- smallprev
		#p100 <- prev*2
		p1 <- smallprev
		p50 <- smallprev
		p100 <- prev*1.5
		
		b1 = runif(23, min=p1, max=p50)
		b2  = runif(50, min=mean(c(p1,p50)), max=mean(c(p50, p100)))
		b3 = runif(25, min=p50, max=p100)
		
		reqdprev <- c(p1, b1, b2, b3, p100)
		
		exdf <- data.frame(x = c(1:100), y = reqdprev)
		pexp <- ggplot(data = exdf, aes(x, y) ) +  xlab("Percentile PRS") + ylab("Expected Prev") + 
		geom_smooth(method = "lm", formula = y ~ poly(x, 3, raw = T)) + geom_point() + theme_bw()
		pexp
		
		polydf <- data.frame()
		for(i in 1:100){
		sub <- dfull[dfull$PGSbin == i,]
		initprev <- sub$Prev[1]
		finalprev <- reqdprev[i]
		sub$AddCase <- 0
		sub$NEWCASE <- sub$CASE
		sub$NEWCASE[sub$NEWCASE == "TBD"] <- "PRSCTRL"
		if(finalprev > initprev){
			sub$AddCase <- (finalprev*sub$N[1]) - sub$CarrCase[1]
			ac <- round((finalprev*sub$N[1]) - sub$CarrCase[1], 0)
			subx1 <- sub[sub$carriers == 1,]
			subx1$NEWCASE <- subx1$CASE
			subx2 <- sub[sub$carriers == 0,]
			subx2$NEWCASE <- "PRSCTRL"
			subx2$NEWCASE[sample(nrow(subx2), ac)] <- "PRSCASE"
			sub <- rbind(subx1, subx2)
		}
		polydf <- rbind(polydf, sub)
		}
		table(polydf$NEWCASE)
		polydf$FINALCASE <- polydf$NEWCASE
		polydf$FINALCASE[polydf$FINALCASE == "PRSCASE"] <- "CASE"
		polydf$FINALCASE[polydf$FINALCASE == "PRSCTRL"] <- "CTRL"
		table(polydf$FINALCASE)
		
		# Prev percentile 
		n=100
		dx <- polydf
		names(dx)[5] <- "invGRS"
		df <- data.frame()
		dfull <- data.frame()
		for(i in 0:(n-1)){
			print(i)
			if(i == 0){
				dx1 <- dx[which(dx$invGRS >= quantile(dx$invGRS, i/n) & dx$invGRS <= quantile(dx$invGRS, (i+1)/n)),]
			}
			else {
				dx1 <- dx[which(dx$invGRS > quantile(dx$invGRS, i/n) & dx$invGRS <= quantile(dx$invGRS, (i+1)/n)),]
			}
			prev3 <- table(dx1$FINALCASE)["CASE"]/nrow(dx1)
			if(is.na(prev3)){prev3 = 0}
			
			# prev perc df
			df2 = data.frame("PGS" = i+1,"meanPGS" = mean(dx1$invGRS) , "Prev3" = prev3) 
			df <- rbind(df, df2)
		}
		pfinal <- df
		p_combined <- ggplot(pfinal, aes(x=PGS, y=Prev3)) + 
		geom_point( size = 2.3) + theme_bw() + xlab("Percentile PRS") + ylab(paste0("Prevalence of disease")) + 
		geom_smooth(method = "lm",formula = y ~ poly(x, 3, raw = T),  size = 1, se = T) +
		ggtitle(paste0("Cases: Carrier + Noncarrier\nTotal n = ", tot))
		
		print("combineddone")
		# Case Ctrls
		names(polydf)[5] <- "invGRS"
		cas <- round(mean(polydf$invGRS[polydf$FINALCASE == "CASE"]), 3)
		ctr <- round(mean(polydf$invGRS[polydf$FINALCASE == "CTRL"]), 3)
		box_full <- ggplot(polydf, aes(y=invGRS, x= FINALCASE, group = as.factor(FINALCASE), fill = as.factor(FINALCASE) )) + geom_boxplot() + 
		ggtitle(paste0("All Individuals n=",tot), 
				subtitle = paste0("Mean PRS CASE = ", cas, "\nMean PRS CTRL = ", ctr))  + 
				xlab("") + ylab("PRS ")  + scale_fill_manual(name = "",values = c("#643843", "grey")) + 
				theme_bw() + theme(legend.position = "None")
		
		# Within case carriers vs noncarriers
		polysub <- polydf[polydf$FINALCASE == "CASE",]
		nc <- round(mean(polysub$invGRS[polysub$carriers == 0]), 3)
		c <- round(mean(polysub$invGRS[polysub$carriers == 1]), 3)
		polysub$CARR[polysub$carriers == 0] <- "NonCarriers"
		polysub$CARR[polysub$carriers == 1] <- "Carriers"
		
		car1 <- table(polysub$carriers)["1"]
		car2 <- table(polysub$carriers)["0"]
		
		t=format(t.test(polysub$invGRS[polysub$carriers == 1], polysub$invGRS[polysub$carriers == 0])$p.value, 3)
		box_cases <- ggplot(polysub, aes(y=invGRS, x= CARR, group = as.factor(CARR), fill = as.factor(CARR) )) + geom_boxplot() + 
		ggtitle("Within Cases", subtitle = paste0("Mean PRS nonCarriers = ", nc, "\nMean PRS Carriers = ", c, "\np-value=",t)) + xlab("") + ylab("PRS of Cases") + 
		scale_fill_manual(values = c("#FC997C", "#396EB0")) + theme_bw() + 
		scale_x_discrete(labels = c(paste0("Carriers\nn=",car1), paste0("NonCarriers\nn=",car2))) + 
		theme(legend.position = "None", 
				axis.text = element_text(size = 12, color = "black"))
		box_cases

		dx <- polysub
		carnum <- table(dx$carriers)["1"]
		noncarnum <- table(dx$carriers)["0"]
		df <- data.frame()
		for(i in 0:(n-1)){
			print(i)
			if(i == 0){
				dx1 <- dx[which(dx$invGRS >= quantile(dx$invGRS, i/n) & dx$invGRS <= quantile(dx$invGRS, (i+1)/n)),]
			}
			else {
				dx1 <- dx[which(dx$invGRS > quantile(dx$invGRS, i/n) & dx$invGRS <= quantile(dx$invGRS, (i+1)/n)),]
			}
			prev_carr <- table(dx1$carriers)["1"]/carnum
			prev_nc <- table(dx1$carriers)["0"]/noncarnum
			if(is.na(prev_carr)){prev_carr = 0}
			if(is.na(prev_nc)){prev_nc = 0}
			
			# prev perc df
			df2 = data.frame("PGS" = i+1, "meanPGS" = mean(dx1$invGRS) , "Prev_carr" = prev_carr, "Prev_NC" = prev_nc) 
			df <- rbind(df, df2)
		}
		
		poly1 <- df[,c("PGS", "meanPGS", "Prev_carr")]
		names(poly1)[3] <- "Prev"
		poly1$G <- "Carriers"
		
		poly2 <- df[,c("PGS", "meanPGS", "Prev_NC")]
		names(poly2)[3] <- "Prev"
		poly2$G <- "NonCarriers"
		
		poly <- rbind(poly1, poly2)
		
		pf <-  ggplot(poly, aes(x=PGS, y=Prev, group = G, color = G)) + 
		geom_point( size = 2.3) + theme_bw() + xlab("Percentile PRS") + ylab(paste0("Proportion of cases")) + 
		geom_smooth(method = "lm",formula = y ~ poly(x, 3, raw = T),  size = 1, se = T)  + 
		scale_color_manual(name = "", values = c("#FC997C", "#396EB0")) + theme(legend.position = "top")
		pf
		
		tot_carriers <- table(polydf$carriers)["1"]
		tot_noncarriers <- table(polydf$carriers)["0"]
		or_carrierrisk <- round((carnum/(tot_carriers-carnum))/(noncarnum/(tot_noncarriers - noncarnum)), 0)
		
		or_PRSrisk <- round(pfinal$Prev3[pfinal$PGS == 100]/pfinal$Prev3[pfinal$PGS == 1],0)
		
		prsdiff <- c-nc
		
		result_df <- data.frame("Total" = tot, "Prevalence" = prev,
								"Total_cases"= nrow(polysub), 
								"CarrierProp"= carrier_perc, 
								"TotalCarriers" = tot_carriers, "Total_NonCarriers"= tot_noncarriers,
								"Penetrance" = car_case_perc, 
								"Carrier_cases" = carnum, "NonCarrier_cases" = noncarnum, 
								"Carrier_case_prop" = carnum/(carnum + noncarnum),
								"OR_Risk_CarriervsNonCarrier" = or_carrierrisk, "OR_PRS_topvsBottom" = or_PRSrisk, 
								"MeanPRSdiff_c_nc" = prsdiff,
								"Prevfinal"= (carnum + noncarnum)/tot)
		resdf <- rbind(resdf, result_df)

		png(paste0("/home/552/sn0095/shared/SarcomaAnalysis/Results/Simulation/SimulationCASES_prev",prev,"_Penet",car_case_perc,"_it",it,"_prs15x.png"), res=350, units = "in", height = 4, width = 8)
		f <- plot_grid(box_cases, pf)
		print({f})
		dev.off()

		png(paste0("/home/552/sn0095/shared/SarcomaAnalysis/Results/Simulation/SimulationALL_prev",prev,"_Penet",car_case_perc,"_it",it,"_prs15x.png"), res=350, units = "in", height = 4, width = 12)
		f1 <- plot_grid(p_carriers, p_combined, box_full, ncol=4)
		print({f1})
		dev.off() 
		print("plots+resdf done")
	}
	resdf$iteration <- it
	write.table(resdf, paste0("/home/552/sn0095/shared/SarcomaAnalysis/Results/Simulation/Simulation_prev",prev,"_Penet",car_case_perc,"_it",it,"fullresults_prs15x.txt"), 
            row.names=F, col.names=T, quote=F, sep = "\t")
	fullres <- rbind(fullres, resdf)
}
# final plot 
r1 <- fullres
write.table(fullres, paste0("/home/552/sn0095/shared/SarcomaAnalysis/Results/Simulation/Simulation_prev",prev,"_Penet",car_case_perc,"fullresults_prs15x.txt"), 
            row.names=F, col.names=T, quote=F, sep = "\t")

f1 <- summarySE(r1, measurevar="MeanPRSdiff_c_nc", groupvars=c("Penetrance"))
f2 <- summarySE(r1, measurevar="OR_Risk_CarriervsNonCarrier", groupvars=c("Penetrance"))
f3 <- summarySE(r1, measurevar="OR_PRS_topvsBottom", groupvars=c("Penetrance"))
f1$OR_carrier <- f2$OR_Risk_CarriervsNonCarrier
f1$PRSfold <- f3$OR_PRS_topvsBottom

write.table(f1, paste0("/home/552/sn0095/shared/SarcomaAnalysis/Results/Simulation/Simulation_prev",prev,"_Penet",car_case_perc,"Summaryresults_prs15x.txt"), 
            row.names=F, col.names=T, quote=F,sep = "\t")

#f1 <- f1[f1$OR_carrier != 0,]
px <- ggplot(f1, aes(x=Penetrance, y=MeanPRSdiff_c_nc)) + geom_point(aes(size=OR_carrier)) + #geom_line() +
  geom_errorbar(aes(ymin=MeanPRSdiff_c_nc-se, ymax=MeanPRSdiff_c_nc+se))  +
  geom_hline(yintercept=0,linetype="dashed", color="black") + 
  xlab("Population Peneterance") + ylab("Mean PRS Carrier - NonCarrier Case") + 
  ggtitle("", subtitle = paste0("n=",tot,", Prev=",prev,", \nCarriers=",carrier_perc,", PRSrisk=2x")) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, size=11, hjust=1),
        axis.text.y = element_text(size=11, hjust=1)) 

png(paste0("/home/552/sn0095/shared/SarcomaAnalysis/Results/Simulation/Simulation_prev",prev,"_Penet",car_case_perc,"Summaryresults_prs15x.png"), res=350, units = "in", height = 4, width = 6)
print({px})
dev.off() 
    
    
