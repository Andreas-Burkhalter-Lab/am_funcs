#ifndef THREADPOOL_HPP
#define THREADPOOL_HPP

#include <pthread.h>
#include <queue>
#include <vector>
using std::queue;
using std::vector;

class Cond;

class Mutex {
   pthread_mutex_t m;
public:
   Mutex(){
      pthread_mutex_init(&m, NULL);
   }
   ~Mutex(){
      pthread_mutex_destroy(&m);
   }
   void lock(){
      pthread_mutex_lock(&m);
   }
   void acquire(){lock();}
   void unlock(){
      pthread_mutex_unlock(&m);
   }
   void release(){unlock();}

   friend class Cond;
};//class, Mutex

class Cond {
   pthread_cond_t c;
public:
   Cond(){
      pthread_cond_init(&c, NULL);
   }
   ~Cond(){
      pthread_cond_destroy(&c);
   }
   void signal(){
      pthread_cond_signal(&c);
   }
   void broadcast(){
      pthread_cond_broadcast(&c);
   }
   void wait(Mutex* mutex){
      pthread_cond_wait(&c, &mutex->m);
   }
};//class, Cond

class Event {
   bool shouldWait;
   Cond wc;
   Mutex mutex;
public:
   Event(){
      shouldWait=true; //initially no signal exists
   }

   //clear/remove any signal
   void reset(){ //todo: rename it
      mutex.lock();
      shouldWait=true;
      mutex.unlock();
   }

   //signal
   void set(){ //todo: rename it
      mutex.lock();
      shouldWait=false;
      wc.broadcast();
      mutex.unlock();
   }

   //wait for the signal
   void wait(){
      mutex.lock();
      while(shouldWait)wc.wait(&mutex);
      mutex.unlock();
   }

   //wait for the signal and reset it
   void waitAndReset(){
      mutex.lock();
      while(shouldWait)wc.wait(&mutex);
      shouldWait=true;
      mutex.unlock();
   }
};//class, Event


typedef void (*Proc)(void* arg);

class ThreadPool;

class Thread {
   Proc p;
   void* arg;
   pthread_t th;
   Event event;
   ThreadPool *pool;

   static void* threadEntry(void* arg);
   friend class ThreadPool;

public:
   Thread(ThreadPool* aPool){
      p=0;
      arg=0;
      pool=aPool;
      pthread_create(&th, NULL, threadEntry, this);
   }//ctor,

   void awaken(){
      event.set();
   }
};//class, Thread


class ThreadPool {
   int N;
   queue<Thread*> que;
   vector<Thread*> allThreads;
   Mutex que_lock;
   Cond que_full_cond;
public:
   bool addTask(Proc p, void* arg){
      if(que.empty()) return false;
      que_lock.acquire();
      if(que.empty()){
         que_lock.release();
         return false;
      }
      Thread* th=que.front();
      que.pop();
      que_lock.release();
      th->p=p;
      th->arg=arg;
      th->awaken();
      return true;
   }
   void addBackIdleThread(Thread* th){
      que_lock.acquire();
      que.push(th);
      if(que.size()==N) que_full_cond.broadcast();
      que_lock.release();
   }
   void wait4all(){
      que_lock.acquire();
      while(que.size()!=N) que_full_cond.wait(&que_lock);
      que_lock.release();
   }

   ThreadPool(int N){
      this->N=N;
      que_lock.acquire();
      for(int i=0; i<N; i++){
	 Thread* th=new Thread(this);
	 que.push(th);
	 allThreads.push_back(th);
      }
      que_lock.release();
   }
};//class, ThreadPool


inline void* Thread::threadEntry(void* arg){
   //todo: maybe we should wait or use mutex on Thread's members?

   Thread* th=(Thread*)arg;
   while(true){
      th->event.waitAndReset();
      th->p(th->arg);
      th->pool->addBackIdleThread(th);
   }//while, true
}//Thread::threadEntry(),

#endif //THREADPOOL_HPP
