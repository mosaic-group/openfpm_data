#!groovy

parallel (

"nyu" : {node ('nyu')
                  {
                    deleteDir()

                    int ntry = 5
                    while (ntry != 0)
                    {
                      try {
                        checkout scm
                        ntry = 0
                      }
                      catch (IOException e)
                      {
                        ntry--
                        sleep(5000)
                      }
                    }

                    stage ('build_nyu')
                    {
                      sh "./build.sh $WORKSPACE $NODE_NAME"
                    }

                    stage ('run_nyu')
                    {
                      sh "cd openfpm_io && ./run.sh $WORKSPACE $NODE_NAME"
                      sh "cd openfpm_io && ./success.sh 2 nyu openfpm_io"
                    }
                  }
                 },


"sb15" : {node ('sbalzarini-mac-15')
                  {
                    deleteDir()
                    env.PATH = "/usr/local/bin:${env.PATH}"

                    int ntry = 5
                    while (ntry != 0)
                    {
                      try {
                        checkout scm
                        ntry = 0
                      }
                      catch (IOException e)
                      {
                        ntry--
                        sleep(5000)
                      }
                    }

                    stage ('build_sb15')
                    {
                      sh "echo $PATH && ./build.sh $WORKSPACE $NODE_NAME"
                    }

                    stage ('run_sb15')
                    {
                      sh "cd openfpm_io && ./run.sh $WORKSPACE $NODE_NAME"
                      sh "cd openfpm_io && ./success.sh 2 sbalzarini-mac-15 openfpm_io"
                    }
                  }
                 },

"gin" : {node ('gin')
                  {
                    deleteDir()

                    int ntry = 5
                    while (ntry != 0)
                    {
                      try {
                        checkout scm
                        ntry = 0
                      }
                      catch (IOException e)
                      {
                        ntry--
                        sleep(5000)
                      }
                    }

                    stage ('build_gin')
                    {
                      sh "echo $PATH && ./build.sh $WORKSPACE $NODE_NAME"
                    }

                    stage ('run_gin')
                    {
                      sh "cd openfpm_io && ./run.sh $WORKSPACE $NODE_NAME"
                      sh "cd openfpm_io && ./success.sh 2 gin openfpm_io"
                    }
                  }
                 }
)

