#!groovy

/* 
   Maintainer name : Mauricio Jimenez 
   Maintainer's email address:  mjjimenezs@eafit.edu.co
 
*/

         /*
         Declare global variable here
         */
         
         // Keep this variable No by default
         DebugMode = "No"
         CustomBranchBuild = DebugMode
        
        
pipeline {


      
        /* 
        Declare the agent where you want to run your build job
        */
        agent any

        /* 
        Tools to use
        
        tools {
        	
        }*/

        /*
        Target environment 
        */
//        parameters {
//            //string defaultValue: 'No', description: '(Not for development team )Choose Yes , if you want to build any custom branch (all branches other than master,develop and feature* are custom branch)' , name: 'CustomBranchBuild', trim: true
//            //choice choices: "No\nYes", description: '(Not for development team )Choose Yes , if you want to build any custom branch (all branches other than master,develop and feature* are custom branch)', name: 'CustomBranchBuild' 
//            //choice choices: "LatestGitCommit\nLatestGitTag", description: 'Choose LatestGitCommit OR LatestGitTag , if you want to build by LatestGitCommit or by LatestGitTag', name: 'CustomBuild'
//
//        }

		parameters {
                   choice choices: "None\nqa\nprod", description: 'Tag selected version to be deployed in the environment', name: 'TagToEnvironment'
                   choice choices: "No\nYes", description: 'Choose "Yes" if selected version is ready to be a release (do tag git) else you may choose No', name: 'IsRelease'
				   choice choices: "develop\nmaster", description: 'Choose "develop" if selected version is ready to be a develop (do tag git) else you may choose master', name: 'branch_name'
        }

        /*
        Print the timestamp in console output
        Discard old build
        Do not allow to run concurrent build 
        */

        options {
                timestamps()
                buildDiscarder logRotator(artifactDaysToKeepStr: '', artifactNumToKeepStr: '', daysToKeepStr: '20', numToKeepStr: '90')
                disableConcurrentBuilds()
        }

        /* 
        Global Environment variable
        */
        environment {
		EMAIL = 'mjjimenezs@eafit.edu.co'
                REPO_URL='gitlab.com/caoba-covid/caoba-eafit/modelo-epidemiologico.git'
        }

        /*
        Stages of pipelines
        */
		
		 stages{
		
                /* 
                This stage is to test block of code , since every time a code is run it takes lot of time to generate results
                and takes lot of time for implementing something new in Jenkinspipeline
                This stage should be uncommented as soon as testing is done.
                */
		         
		        stage('Test Certain block') {
		
		               steps {
		                               
		                     sh '''
		                     echo "Blue Ocean path is : ${RUN_DISPLAY_URL}"
		                     ''' 
		                                          
		               }
		        }
                /* 
                Code checkout stage
                */
                stage('Checkout') {
                        when {
							beforeAgent true
							branch 'develop'							
                        	 // check the or !   
                             expression { ( branch_name =~ /^feature.*/ ||  branch_name =~ /^hotfix.*/ || branch_name =~ /^bugfix.*/ || branch_name == 'master' ||  (branch_name == 'develop' ||  CustomBranchBuild == 'Yes' ) ) }
                        }
                        steps {
                                checkout scm
                                                       
                                script {
                                        /*
		                                Variable to control the flow of job
		                                */
		                                env.IS_AUTOMATED_COMMIT="NO"
				                        env.GO_AHEAD="YES"
		                                
                                		lastCommitMessage=sh(script: "git log --format=%B -n 1", returnStdout: true).trim()
                                		print("Last commit message" + lastCommitMessage)
                                        if (lastCommitMessage.startsWith("Commit by JenkinsServer")) {
                                                 echo "This is an increment pom push done by JenkinsServer, hence omitting all build stages subsequently"
                                                 env.IS_AUTOMATED_COMMIT =  "YES"
                                                 print("auto commit " + IS_AUTOMATED_COMMIT)
                                                 if( IS_AUTOMATED_COMMIT ==  "YES")
                                                    {       
    													GO_AHEAD="NO"
											    	}
                                        } else {
                                           echo "Commit not pushed by Jenkinsserver , hence execution will be done for all stages"
                                        }

                                      
                                }             
                        }
                }
				
                stage('Upload trusted zone to Datalake S3') {                      
                        when {             
                              expression { (branch_name == 'master' ||  (branch_name == 'develop' ||  CustomBranchBuild == 'Yes' ) ) }
                        }

                        steps {
                               sh 'echo upload....'
                               dir('/var/lib/jenkins/workspace/modelo-epidemiologico/data') {
                                        pwd(); //Log current directory

                                        withAWS(region:'us-east-2', credentials:'e915e1b8-ce3c-4cc1-bfe3-033db85278bc') {

                                                //def identity=awsIdentity();//Log AWS credentials

                                                // Upload files from working directory 'dist' in your project workspace
                                                s3Upload(bucket:"trusted-epidemiologico", workingDir:'trusted', includePathPattern:'**/*');
                                        }

                                };
                                sh 'echo done upload'
                        }
                }

                stage('Upload refined zone to Datalake S3') {                      
                        when {             
                              expression { (branch_name == 'master' ||  (branch_name == 'develop' ||  CustomBranchBuild == 'Yes' ) ) }
                        }

                        steps {
                              sh 'echo upload....'
                               dir('/var/lib/jenkins/workspace/modelo-epidemiologico/data') {
                                        pwd(); //Log current directory

                                        withAWS(region:'us-east-2', credentials:'e915e1b8-ce3c-4cc1-bfe3-033db85278bc') {

                                                //def identity=awsIdentity();//Log AWS credentials

                                                // Upload files from working directory 'dist' in your project workspace
                                                s3Upload(bucket:"refined-epidemiologico", workingDir:'refined', includePathPattern:'**/*');
                                        }

                                };
                                sh 'echo done upload'
                        }
                }
        }


        /* 
        Post stage
        */

        post {

                success {
                        /*
                        Post success notification
                        */
                        script {
                               if(env.SHOULD_PUSH=="YES"){
                                    //try {

							  	    //withCredentials([usernamePassword(credentialsId: 'mjjimenezs', passwordVariable: 'GIT_PASSWORD', usernameVariable: 'GIT_USERNAME')]) {
									//sh('echo " The current branch from script is   -> $CURRENT_BRANCH" ')
									//sh('git push https://${GIT_USERNAME}:${GIT_PASSWORD}@$REPO_URL ')
								  //sh('git push --tags https://${GIT_USERNAME}:${GIT_PASSWORD}@$REPO_URL ')		  
                                    //}
                               //}
                               //catch(Exception ex) {
								//              println(ex);
							   //}
                                   
                               }
								echo "Success deploy"
							}         
                 emailext (
                        	body: "From ${env.JOB_NAME} : Job ${env.BUILD_ID} build success , Normal Url : ${env.BUILD_URL} , BlueOceanUrl : ${env.RUN_DISPLAY_URL} , SonarResults : ${env.SonarResultUrl}",
                                replyTo: '${EMAIL}',
                                subject: '${DEFAULT_SUBJECT}',
                                to: emailextrecipients(
                                	[
                                         	[
                                                	$class: 'CulpritsRecipientProvider'
                                                ],
                                                [
                                                        $class: 'RequesterRecipientProvider'
                                                ]
                                	]
                        	)
                	)
                }
                
                /* 
                If build fails mail to focal point
                */
                failure {
                        /*
                        Post failure notification
                        */
                 emailext (
                        	body: "From ${env.JOB_NAME} : Job ${env.BUILD_ID} build failure!! , Normal Url : ${env.BUILD_URL} , BlueOceanUrl : ${env.RUN_DISPLAY_URL} , SonarResults : ${env.SonarResultUrl}",
                                replyTo: '${EMAIL}',
                                subject: '${DEFAULT_SUBJECT}',
                                to: emailextrecipients(
                                	[
                                         	[
                                                	$class: 'CulpritsRecipientProvider'
                                                ],
                                                [
                                                        $class: 'RequesterRecipientProvider'
                                                ]
                                	]
                        	)
                	)
                }
        }

}
