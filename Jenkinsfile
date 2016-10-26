#!groovy

parallel (

"nyu_VALGRIND" : {node ('nyu')
                  {
                    deleteDir()
                    checkout scm
                    stage ('build_nyu_val')
                    {
                      sh "./build.sh $WORKSPACE $NODE_NAME VALGRIND"
                    }

                    stage ('run_nyu_val')
                    {
                      sh "cd openfpm_data && ./run.sh $WORKSPACE $NODE_NAME VALGRIND"
                    }
                  }
                 },


"nyu_NO" : {node ('nyu')
                  {
                    deleteDir()
                    checkout scm
                    stage ('build_nyu_nor')
                    {
                      sh "./build.sh $WORKSPACE $NODE_NAME NO"
                    }

                    stage ('run_nyu_nor')
                    {
                      sh "cd openfpm_data && ./run.sh $WORKSPACE $NODE_NAME NO"
                    }
                  }
                 },

"nyu_SE" : {node ('nyu')
                  {
                    deleteDir()
                    checkout scm
                    stage ('build_nyu_se')
                    {
                      sh "./build.sh $WORKSPACE $NODE_NAME SE"
                    }

                    stage ('run_nyu_se')
                    {
                      sh "cd openfpm_data && ./run.sh $WORKSPACE $NODE_NAME SE"
                    }
                  }
                 },

"sb15_VALGRIND" : {node ('sbalzarini-mac-15')
                  {
                    deleteDir()
                    env.PATH = "/usr/local/bin:${env.PATH}"
                    checkout scm
                    stage ('build_sb15_val')
                    {
                      sh "echo $PATH && ./build.sh $WORKSPACE $NODE_NAME VALGRIND"
                    }

                    stage ('run_sb15_val')
                    {
                      sh "cd openfpm_data && ./run.sh $WORKSPACE $NODE_NAME VALGRIND"
                    }
                  }
                 },


"sb15_NO" : {node ('sbalzarini-mac-15')
                  {
                    deleteDir()
                    env.PATH = "/usr/local/bin:${env.PATH}"
                    checkout scm
                    stage ('build_sb15_nor')
                    {
                      sh "./build.sh $WORKSPACE $NODE_NAME NO"
                    }

                    stage ('run_sb15_nor')
                    {
                      sh "cd openfpm_data && ./run.sh $WORKSPACE $NODE_NAME NO"
                    }
                  }
                 },

"sb15_SE" : {node ('sbalzarini-mac-15')
                  {
                    deleteDir()
                    env.PATH = "/usr/local/bin:${env.PATH}"
                    checkout scm
                    stage ('build_sb15_se')
                    {
                      sh "./build.sh $WORKSPACE $NODE_NAME SE"
                    }

                    stage ('run_sb15_se')
                    {
                      sh "cd openfpm_data && ./run.sh $WORKSPACE $NODE_NAME SE"
                    }
                  }
                 }



)

